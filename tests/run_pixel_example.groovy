// CLI with Fiji libraries: run_pixel_example.groovy RELEASE NEW_OUTPUT_FOLDER
import groovy.json.JsonSlurper
import java.nio.file.Files
import java.nio.file.StandardCopyOption
try {
    File root=new File(args[0]),out=new File(args[1])
    if(out.exists())throw new IllegalStateException('Output already exists')
    File input=new File(root,'examples/input')
    input.eachFileRecurse{f->
        if(f.isFile()){
            File dest=new File(out,input.toPath().relativize(f.toPath()).toString())
            dest.parentFile.mkdirs();Files.copy(f.toPath(),dest.toPath())
        }
    }
    File exact=new GroovyShell(this.class.classLoader,new Binding([batchDirectory:out,minPunctumAreaUm2:0.02d,maxPunctumAreaUm2:400d])).evaluate(new File(root,'scripts/exact_pixels.groovy'))
    def parser=new JsonSlurper(),actual=parser.parse(new File(exact,'Geometry_QC_results.json')),expected=parser.parse(new File(root,'examples/expected/exact_regions.json'))
    assert actual.regions.size()==expected.Regions
    actual.regions.eachWithIndex{r,i->
        def e=expected.regions[i]
        assert r.Unique_ID==e.Unique_ID
        ['Area_pixels','Area_um2','Integrated_density_raw','Mean_raw_intensity','Max_raw_intensity'].each{k->assert Math.abs((r[k] as double)-(e[k] as double))<1e-9}
    }
    def label=ij.IJ.openImage(new File(exact,'Labels/0017_hole_aware.tif').path)
    def target=ij.IJ.openImage(new File(root,'examples/expected/0017_hole_aware.tif').path)
    assert java.util.Arrays.equals((short[])label.processor.pixels,(short[])target.processor.pixels)
    label.close();target.close()
    println('PASS: 58 exact regions, 2 corrected geometries, identical label pixels.');System.exit(0)
}catch(Throwable e){e.printStackTrace();System.exit(1)}
