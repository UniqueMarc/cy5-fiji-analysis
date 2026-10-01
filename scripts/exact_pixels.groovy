import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import ij.IJ
import ij.ImagePlus
import ij.Prefs
import ij.gui.Roi
import ij.io.RoiDecoder
import ij.io.RoiEncoder
import ij.io.FileSaver
import ij.measure.Measurements
import ij.measure.ResultsTable
import ij.plugin.filter.ParticleAnalyzer
import ij.process.ImageStatistics
import ij.process.ImageProcessor
import ij.process.ShortProcessor
import java.util.zip.ZipFile
import java.util.zip.ZipOutputStream
import java.util.zip.ZipEntry

// Same ImageJ particle scan, thresholds, connectivity and area limits as the
// original analysis. The only changed option is COMPOSITE_ROIS. This preserves
// internal holes and excludes separately detected islands from outer ROIs.
@CompileStatic
class CollectCompositeParticles extends ParticleAnalyzer {
 List<Roi> shapes=new ArrayList<Roi>()
 List<Integer> componentPixels=new ArrayList<Integer>()
 List<int[]> seeds=new ArrayList<int[]>()
 CollectCompositeParticles(double lo,double hi){super(ParticleAnalyzer.SHOW_NONE | ParticleAnalyzer.COMPOSITE_ROIS,Measurements.AREA,new ResultsTable(),lo,hi,0d,1d)}
 @Override protected void saveResults(ImageStatistics st,Roi roi){
  shapes.add((Roi)roi.clone());componentPixels.add(st.pixelCount);seeds.add([st.xstart,st.ystart] as int[])
 }
}
@CompileStatic
class CheckGeometry {
 static Map measure(Roi oldRoi,Roi newRoi,ImageProcessor raw,ImageProcessor binary,int[] oldCounts,short[] correctedLabels,int label,int seedX,int seedY) {
  int w=raw.getWidth();java.awt.Rectangle b=oldRoi.getBounds()
  int bx=(int)b.getX(),by=(int)b.getY(),bw=(int)b.getWidth(),bh=(int)b.getHeight()
  long oldSum=0,newSum=0,sx=0,sy=0;int oldN=0,newN=0,oldBackground=0,removedForeground=0,maxVal=0
  java.awt.Rectangle nb=newRoi.getBounds()
  if(!b.equals(nb))throw new IllegalStateException('Outer bounds changed for label '+label)
  for(int y=by;y<by+bh;y++)for(int x=bx;x<bx+bw;x++){
   boolean a=oldRoi.contains(x,y);int p=y*w+x,v=raw.get(p)
   if(a){oldN++;oldSum+=v;oldCounts[p]++;if(binary.get(p)==0)oldBackground++}
  }
  // Work directly with the archived binary component, not a vector ROI round
  // trip. Validate the flood-filled pixel count against ImageJ's internal
  // ParticleAnalyzer count. ShapeRoi.contains may differ at thin contours.
  int h=raw.getHeight(),seed=seedY*w+seedX
  if(binary.get(seed)==0||correctedLabels[seed]!=0)throw new IllegalStateException('Invalid/duplicate component seed '+label)
  int[] queue=new int[w*h];int head=0,tail=0;queue[tail++]=seed;correctedLabels[seed]=(short)label
  while(head<tail){int p=queue[head++],x=p%w,y=(int)(p/w),v=raw.get(p)
   if(!oldRoi.contains(x,y))throw new IllegalStateException('Component extends beyond original outer polygon '+label)
   newN++;newSum+=v;sx+=x;sy+=y;maxVal=Math.max(maxVal,v)
   for(int dy=-1;dy<=1;dy++)for(int dx=-1;dx<=1;dx++){
    int xx=x+dx,yy=y+dy;if(xx<0||yy<0||xx>=w||yy>=h)continue
    int q=yy*w+xx;if(binary.get(q)==0)continue
    if(correctedLabels[q]!=0&&correctedLabels[q]!=label)throw new IllegalStateException('Duplicate components '+label)
    if(correctedLabels[q]==0){correctedLabels[q]=(short)label;queue[tail++]=q}
   }
  }
  for(int y=by;y<by+bh;y++)for(int x=bx;x<bx+bw;x++)if(oldRoi.contains(x,y)&&binary.get(x,y)>0&&correctedLabels[y*w+x]!=label)removedForeground++
  return [Original_pixels:oldN,Original_sum:oldSum,Corrected_pixels:newN,Corrected_sum:newSum,Original_background_pixels:oldBackground,Removed_foreground_other_components:removedForeground,Max_raw_intensity:maxVal,X_pixel:((double)sx)/newN,Y_pixel:((double)sy)/newN]
 }
}
// batchDirectory is a stage-1 output; never overwrite a previous correction.
File batch = batchDirectory
File out = new File(batch,'Exact_pixel_measurements')
if(out.exists()) throw new IllegalStateException('Output exists; will not overwrite: '+out)
out.mkdirs(); new File(out,'Labels').mkdirs(); new File(out,'Corrected_ROIs').mkdirs()
def parser=new JsonSlurper(),d=parser.parse(new File(batch,'All_field_results.json'))
def regions=d.regions.groupBy{it.Image_Index}
List rows=[],images=[]
def selected=d.images
Prefs.blackBackground=true
selected.each{im->
 int idx=im.Image_Index as int;def ident=[Current_group:im.Group,Current_ND2_path:im.Source_path]
 File folder=new File(batch,'640_Cy5_ZMAX/'+im.Source_folder)
 if(!folder.isDirectory())throw new IllegalStateException('Missing source folder')
 File rawFile=new File(folder,im.Image+'_640_Cy5_ZMAX_RAW16.tif')
 File maskFile=new File(batch,'Cy5_Puncta_Masks/'+folder.name+'/'+im.Image+'_Cy5_ZMAX_puncta_mask.tif')
 File oldRoiFile=new File(batch,'ROI_Files/'+folder.name+'/'+im.Image+'_ALL_Cy5_ROIs.zip')
 def raw=IJ.openImage(rawFile.path),mask=IJ.openImage(maskFile.path)
 if(raw==null||mask==null||raw.width!=mask.width||raw.height!=mask.height)throw new IllegalStateException('Missing/incorrect raw or mask image '+idx)
 double pixelArea=(im.Pixel_width_um as double)*(im.Pixel_height_um as double)
 def analyzer=new CollectCompositeParticles(minPunctumAreaUm2/pixelArea,maxPunctumAreaUm2/pixelArea)
 def work=mask.duplicate();work.processor.setThreshold(255,255,ImageProcessor.NO_LUT_UPDATE)
 if(!analyzer.analyze(work))throw new IllegalStateException('Particle analysis failed '+idx)
 List<Roi> oldRois=[];ZipFile zip=oldRoiFile.isFile()?new ZipFile(oldRoiFile):null
 try{if(zip!=null)zip.entries().each{e->if(e.name.endsWith('.roi')){def stream=zip.getInputStream(e);try{oldRois.add(new RoiDecoder(stream.bytes,e.name).getRoi())}finally{stream.close()}}}}finally{if(zip!=null)zip.close()}
 def originalRows=regions[im.Image_Index] ?: []
 if(analyzer.shapes.size()>65535)throw new IllegalStateException("Too many labels for UINT16 label image")
 if(oldRois.size()!=originalRows.size()||analyzer.shapes.size()!=oldRois.size())throw new IllegalStateException('Particle population changed '+idx+' old='+oldRois.size()+' new='+analyzer.shapes.size())
 int[] oldCounts=new int[raw.width*raw.height];short[] labels=new short[oldCounts.length]
 File roiPath=new File(out,'Corrected_ROIs/'+String.format('%04d',idx)+'_hole_aware.zip')
 ZipOutputStream zo=new ZipOutputStream(new FileOutputStream(roiPath))
 long totalOldSum=0,totalNewSum=0,totalOldPixels=0,totalNewPixels=0;int changed=0
 try{
  for(int j=0;j<oldRois.size();j++){
   def orig=originalRows[j];Roi oldR=oldRois[j],newR=analyzer.shapes[j]
   def seed=analyzer.seeds[j]
   def m=CheckGeometry.measure(oldR,newR,raw.processor,mask.processor,oldCounts,labels,j+1,seed[0],seed[1])
   if(m.Original_pixels!=(orig.Area_pixels as int)||Math.abs((m.Original_sum as double)-(orig.Integrated_density_raw as double))>1e-5)throw new IllegalStateException('Original measurement/identity mismatch '+orig.Unique_ID)
   if(m.Corrected_pixels!=analyzer.componentPixels[j])throw new IllegalStateException('Composite ROI differs from ParticleAnalyzer component '+orig.Unique_ID)
   boolean c=m.Original_pixels!=m.Corrected_pixels;if(c)changed++
   rows.add([Group:ident.Current_group,Image_Index:idx,Region_ID:orig.Region_ID,Unique_ID:orig.Unique_ID,Original_area_um2:orig.Area_um2,Original_integrated_density_raw:orig.Integrated_density_raw,Original_pixels:orig.Area_pixels,Area_um2:(m.Corrected_pixels as double)*pixelArea,Area_pixels:m.Corrected_pixels,Integrated_density_raw:m.Corrected_sum,Mean_raw_intensity:(m.Corrected_sum as double)/(m.Corrected_pixels as double),Max_raw_intensity:m.Max_raw_intensity,Removed_background_pixels:m.Original_background_pixels,Removed_other_component_pixels:m.Removed_foreground_other_components,Changed:c,X_pixel:m.X_pixel,Y_pixel:m.Y_pixel,Current_raw_path:ident.Current_ND2_path])
   totalOldSum+=m.Original_sum;totalNewSum+=m.Corrected_sum;totalOldPixels+=m.Original_pixels;totalNewPixels+=m.Corrected_pixels
   newR.setName(orig.Unique_ID as String);zo.putNextEntry(new ZipEntry(orig.Unique_ID+'.roi'));new RoiEncoder(zo).write(newR);zo.closeEntry()
  }
 }finally{zo.close()}
 int overlapPixels=0,extraPixelMeasurements=0;long duplicatedIntensity=0;int foreground=0;long uniqueSum=0
 for(int p=0;p<labels.length;p++){
  if(oldCounts[p]>1){overlapPixels++;extraPixelMeasurements+=oldCounts[p]-1;duplicatedIntensity+=(long)(oldCounts[p]-1)*raw.processor.get(p)}
  if(labels[p]!=0){foreground++;uniqueSum+=raw.processor.get(p)}
 }
 if(foreground!=totalNewPixels||uniqueSum!=totalNewSum)throw new IllegalStateException('Whole-field uniqueness validation failed '+idx)
 File labelPath=new File(out,'Labels/'+String.format('%04d',idx)+'_hole_aware.tif')
 def outputImage=new ImagePlus('Hole-aware region IDs',new ShortProcessor(raw.width,raw.height,labels,null));outputImage.calibration=raw.calibration.copy()
 if(!new FileSaver(outputImage).saveAsTiff(labelPath.path))throw new IllegalStateException('Save failed '+idx)
 images.add([Group:ident.Current_group,Image_Index:idx,Regions:originalRows.size(),Changed_regions:changed,Nuclei:im.Nuclei,Crops:im.Crops,Original_area_pixels:totalOldPixels,Corrected_area_pixels:totalNewPixels,Original_intensity:totalOldSum,Corrected_intensity:totalNewSum,Original_overlap_pixels:overlapPixels,Original_extra_pixel_measurements:extraPixelMeasurements,Original_duplicated_intensity:duplicatedIntensity,Corrected_overlap_pixels:0,Current_raw_path:ident.Current_ND2_path,Source_raw_ZMAX:rawFile.path,Source_mask:maskFile.path,Corrected_ROI_path:roiPath.path,Corrected_label_path:labelPath.path,Threshold_raw:im.Threshold_raw,Pixel_area_um2:pixelArea])
  raw.close();mask.close();work.close();outputImage.close()
  println('GEOMETRY '+images.size()+'/'+selected.size()+' image='+idx+' regions='+originalRows.size()+' changed='+changed+' old_overlap_px='+overlapPixels)
}
def result=[Parameters:[Only_algorithm_change:'Remeasure exact 8-connected archived-mask components using ParticleAnalyzer seeds and validated internal pixel counts. Detection, acceptance, thresholds and region IDs unchanged.',Vector_ROIs:'COMPOSITE_ROIS exports are visual aids only; canonical measurement geometry is the label TIFF. Vector-to-raster conversion can differ at boundary pixels.',Minimum_area_um2:minPunctumAreaUm2,Maximum_area_um2:maxPunctumAreaUm2,Population:'All original accepted particles, 1:1 stable IDs; raw ZMAX unsmoothed measurements on exact connected foreground.',Correction:'Exclude internal background holes and separately detected internal components from each outer particle ROI. Do not just intersect a filled outer polygon with the binary mask.',Original_data_preserved:true],regions:rows,images:images]
new File(out,'Geometry_QC_results.json').setText(JsonOutput.toJson(result),'UTF-8')
println('FINISHED '+out.path+' regions='+rows.size()+' changed='+rows.count{it.Changed})
return out
