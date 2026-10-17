import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import ij.IJ
import ij.ImagePlus
import ij.gui.Roi
import ij.io.RoiDecoder
import ij.io.FileSaver
import ij.process.ShortProcessor
import java.util.zip.ZipFile

// Preserve original contour centroids for ownership; exact component pixel
// centroids describe geometry and must not replace assignment centroids.
@CompileStatic
class CellTerritories {
    static short[] owners(List<Roi> nuclei, List<Roi> expanded, int w, int h) {
        short[] labels=new short[w*h]
        float[] best=new float[w*h]
        Arrays.fill(best,Float.POSITIVE_INFINITY)
        for (int j=0;j<expanded.size();j++) {
            Roi roi=expanded[j]
            double[] center=nuclei[j].getContourCentroid()
            def b=roi.getBounds()
            for (int y=Math.max(0,(int)b.y);y<Math.min(h,b.y+b.height);y++) {
                for (int x=Math.max(0,(int)b.x);x<Math.min(w,b.x+b.width);x++) {
                    if (!roi.contains(x,y)) continue
                    int p=y*w+x
                    double dx=x-center[0],dy=y-center[1]
                    float distance=(float)(dx*dx+dy*dy)
                    // Equal distances retain the earlier nuclear ROI ID.
                    if (distance<best[p]) {best[p]=distance;labels[p]=(short)(j+1)}
                }
            }
        }
        return labels
    }
}

List<Roi> readCellRois(File file, int expected) {
    if (expected==0 && !file.isFile()) return []
    if (!file.isFile()) throw new IllegalStateException('Missing ROI file: '+file)
    List<Roi> rois=[]
    ZipFile zip=new ZipFile(file)
    try {
        zip.entries().each {entry ->
            if (entry.name.endsWith('.roi')) {
                def stream=zip.getInputStream(entry)
                try {rois.add(new RoiDecoder(stream.bytes,entry.name).getRoi())}
                finally {stream.close()}
            }
        }
    } finally {zip.close()}
    return rois
}

File batch=batchDirectory
File exactFolder=new File(batch,'Exact_pixel_measurements')
File output=new File(batch,'Per_cell_measurements')
if (output.exists()) throw new IllegalStateException('Output exists; will not overwrite: '+output)
def parser=new JsonSlurper()
def source=parser.parse(new File(batch,'All_field_results.json'))
def exact=parser.parse(new File(exactFolder,'Geometry_QC_results.json'))
def byImage=exact.regions.groupBy {it.Image_Index as int}
def imageMeta=exact.images.collectEntries {[(it.Image_Index as int):it]}
if (imageMeta.size()!=exact.images.size() || imageMeta.keySet()!=source.images.collect {it.Image_Index as int}.toSet())
    throw new IllegalStateException('Image IDs do not match between pipeline stages')
output.mkdirs()
new File(output,'Territory_labels').mkdirs()
List images=[],cells=[],puncta=[]
source.images.sort {it.Image_Index as int}.each {im ->
    int idx=im.Image_Index as int
    def measured=imageMeta[idx]
    def rows=(byImage[idx] ?: []).sort {it.Region_ID as int}
    File roiRoot=new File(batch,'ROI_Files/'+im.Source_folder)
    List<Roi> nuclei=readCellRois(new File(roiRoot,im.Image+'_nucleus_ROIs.zip'),im.Nuclei as int)
    List<Roi> expanded=readCellRois(new File(roiRoot,im.Image+'_perinuclear_analysis_ROIs.zip'),im.Nuclei as int)
    List<Roi> spots=readCellRois(new File(roiRoot,im.Image+'_ALL_Cy5_ROIs.zip'),rows.size())
    if (nuclei.size()!=expanded.size() || nuclei.size()!=(im.Nuclei as int) || spots.size()!=rows.size())
        throw new IllegalStateException('ROI population mismatch: '+idx)
    if (nuclei.size()>65535) throw new IllegalStateException('Too many nuclear IDs for UINT16 labels')
    def labelImage=IJ.openImage(new File(exactFolder,'Labels/'+String.format('%04d_hole_aware.tif',idx)).path)
    if (labelImage==null) throw new IllegalStateException('Missing exact component labels: '+idx)
    int w=labelImage.width,h=labelImage.height
    short[] owners=CellTerritories.owners(nuclei,expanded,w,h)
    long[] sizes=new long[nuclei.size()+1]
    for (int p=0;p<owners.length;p++) sizes[owners[p]&65535]++
    int assigned=0
    rows.eachWithIndex {row,j ->
        if ((row.Region_ID as int)!=j+1) throw new IllegalStateException('Unstable component IDs: '+idx)
        double[] center=spots[j].getContourCentroid()
        int x=(int)Math.round(center[0]),y=(int)Math.round(center[1])
        int owner=(x>=0 && x<w && y>=0 && y<h) ? (owners[y*w+x]&65535) : 0
        if (owner>0) assigned++
        puncta.add([Image_Index:idx,Region_ID:row.Region_ID,Unique_ID:row.Unique_ID,Cell_ID:owner,
            Area_pixels:row.Area_pixels,Area_um2:row.Area_um2,Integrated_density_raw:row.Integrated_density_raw,
            Mean_raw_intensity:row.Mean_raw_intensity,Assignment_X_pixel:center[0],Assignment_Y_pixel:center[1],
            Component_X_pixel:row.X_pixel,Component_Y_pixel:row.Y_pixel])
    }
    double px=im.Pixel_width_um as double,py=im.Pixel_height_um as double
    nuclei.eachWithIndex {roi,j ->
        double[] center=roi.getContourCentroid()
        cells.add([Image_Index:idx,Cell_ID:j+1,Center_X_pixel:center[0],Center_Y_pixel:center[1],
            Nucleus_area_um2:roi.getStatistics().area*px*py,Territory_pixels:sizes[j+1],Territory_area_um2:sizes[j+1]*px*py])
    }
    def territory=new ImagePlus('Nucleus-associated territory IDs',new ShortProcessor(w,h,owners,null))
    territory.calibration=labelImage.calibration.copy()
    if (!new FileSaver(territory).saveAsTiff(new File(output,'Territory_labels/'+String.format('%04d.tif',idx)).path))
        throw new IllegalStateException('Cannot save territory labels: '+idx)
    images.add([Image_Index:idx,Group:measured.Group,Raw_filename:new File(im.Source_path as String).name,
        Width:w,Height:h,Pixel_width_um:px,Pixel_height_um:py,Nuclear_regions:nuclei.size(),
        All_puncta:rows.size(),Assigned_puncta:assigned,Unassigned_puncta:rows.size()-assigned])
    labelImage.close();territory.close()
    println('PER_CELL image='+idx+' cells='+nuclei.size()+' assigned='+assigned+'/'+rows.size())
}
def result=[Schema_version:2,Parameters:[Territory:'Nuclear ROI expanded outward by 20 micrometres, including nuclear interiors; overlap goes to nearest nuclear contour centroid.',
    Assignment:'Original punctum contour centroid rounded with Math.round; entire punctum assigned once; ID 0 means unassigned.',
    Measurement:'Exact connected-component area and unscaled raw ZMAX integrated fluorescence.',
    Zero_puncta:'Counts and sums are zero; punctum means are null.'],images:images,cells:cells,puncta:puncta]
new File(output,'per_cell_measurements.json').setText(JsonOutput.toJson(result),'UTF-8')
return output
