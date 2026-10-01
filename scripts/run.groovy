import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import groovy.lang.Binding
import groovy.lang.GroovyShell

// GUI supplies the same variables that the CLI fills here.
boolean cli = binding.hasVariable('args') && args.length > 0
if (cli) {
    if(args.length < 5) throw new IllegalArgumentException('Usage: run.groovy RELEASE INPUT_ROOT OUTPUT_PARENT CONFIG_JSON manifest|discover [stable IDs, comma-separated]')
    binding.setVariable('releaseFolder',new File(args[0]))
    binding.setVariable('inputFolder',new File(args[1]))
    binding.setVariable('outputFolder',new File(args[2]))
    binding.setVariable('configFile',new File(args[3]))
    if(!(args[4] in ['manifest','discover'])) throw new IllegalArgumentException('Mode must be manifest or discover')
    binding.setVariable('useFinalManifest',args[4]=='manifest')
}
def perform = {
    def parser=new JsonSlurper()
    Map config=parser.parse(configFile)
    Map p=new LinkedHashMap(config.parameters)
    ['hoechstChannel','cy5Channel','singleCellOutputPixels'].each{p[it]=((Number)p[it]).intValue()}
    p.each{k,v->if(v instanceof Number && !(v instanceof Integer))p[k]=((Number)v).doubleValue()}
    if(p.splitConnectedCy5)throw new IllegalArgumentException('This release implements connected components without peak splitting; splitConnectedCy5 must be false.')
    if(!(p.fixedCy5DisplayMaxRaw > p.fixedCy5DisplayMinRaw && p.fixedCy5DisplayMinRaw>=0))throw new IllegalArgumentException('Specify valid fixed Cy5 display bounds.')
    if(p.cy5DetectionSigmaPixels<0 || p.minPunctumAreaUm2<0 || p.maxPunctumAreaUm2<=p.minPunctumAreaUm2)throw new IllegalArgumentException('Invalid detection parameters')
    assert inputFolder.isDirectory():'Input folder not found'
    assert releaseFolder.isDirectory():'Release folder not found'
    if(!outputFolder.exists())assert outputFolder.mkdirs():'Cannot create output parent'
    List<String> groups=config.group_order.collect{it.toString()}
    if(groups.toSet().size()!=groups.size() || groups.any{!(it==~/[A-Za-z0-9._-]+/)})throw new IllegalArgumentException('Group names must be unique ASCII folder names')
    Map<String,String> labels=[:]
    groups.eachWithIndex{g,i->labels[g]=String.format('%02d_%s',i+1,g)}
    List<Map> entries=[]
    Set selectedIds=cli && args.length>5 ? args[5].split(',').collect{it as int}.toSet() : null
    if(useFinalManifest){
        List records=parser.parse(new File(releaseFolder,'data_manifest.json'))
        records.findAll{it.Included==1 && (selectedIds==null || selectedIds.contains(it.Image_Index as int))}.each{im->
            assert im.Group in groups:'Manifest group missing from config'
            File f=new File(inputFolder,im.Raw_relative_path)
            assert f.canonicalPath.startsWith(inputFolder.canonicalPath+File.separator):'Input outside root'
            assert f.isFile():'Missing manifest input: '+f
            entries.add([file:f,rootLabel:labels[im.Group],rootPath:new File(inputFolder,im.Group).path,relativePath:f.name,stableId:im.Image_Index as int,group:im.Group])
        }
        if(selectedIds!=null && entries.collect{it.stableId}.toSet()!=selectedIds)throw new IllegalArgumentException('Some selected IDs are not in the valid-data manifest')
    }else{
        if(selectedIds!=null)throw new IllegalArgumentException('Stable IDs require manifest mode')
        groups.each{g->
            File dir=new File(inputFolder,g)
            if(dir.isDirectory()){
                List<File> files=[]
                def collectFile={f->if(f.isFile() && f.name.toLowerCase().endsWith('.nd2'))files.add(f)}
                if(p.includeSubfolders)dir.eachFileRecurse(collectFile)
                else dir.eachFile(collectFile)
                files.sort{it.canonicalPath}.each{f->
                    entries.add([file:f,rootLabel:labels[g],rootPath:dir.path,relativePath:dir.toPath().relativize(f.toPath()).toString(),stableId:entries.size()+1,group:g])
                }
            }
        }
    }
    if(entries.isEmpty())throw new IllegalArgumentException('No input files selected')
    assert entries.collect{it.file.canonicalPath}.toSet().size()==entries.size():'Duplicate file'
    p.outputFolder=outputFolder
    p.selectedInputFolders=entries.collect{new File(it.rootPath)}.unique()
    p.suppliedInputEntries=entries.sort{it.stableId}
    loci.common.DebugTools.setRootLevel('WARN')
    Binding context=new Binding(p)
    File batch=new GroovyShell(this.class.classLoader,context).evaluate(new File(releaseFolder,'scripts/pipeline_core.groovy'))
    if(batch==null || !new File(batch,'All_field_results.json').isFile())throw new IllegalStateException('Stage 1 did not complete')
    new File(batch,'Run_configuration.json').setText(JsonOutput.prettyPrint(JsonOutput.toJson(config)),'UTF-8')
    new File(batch,'Run_input_map.json').setText(JsonOutput.prettyPrint(JsonOutput.toJson(entries.collect{[Image_Index:it.stableId,Group:it.group,Source_path:it.file.path]})),'UTF-8')
    def environment=[ImageJ:ij.IJ.getVersion(),BioFormats:loci.formats.FormatTools.VERSION,Groovy:GroovySystem.version,Java:System.getProperty('java.version'),OS:System.getProperty('os.name')]
    new File(batch,'Runtime_versions.json').setText(JsonOutput.prettyPrint(JsonOutput.toJson(environment)),'UTF-8')
    new GroovyShell(this.class.classLoader,new Binding([batchDirectory:batch,minPunctumAreaUm2:p.minPunctumAreaUm2,maxPunctumAreaUm2:p.maxPunctumAreaUm2])).evaluate(new File(releaseFolder,'scripts/exact_pixels.groovy'))
    println('RELEASE_OUTPUT '+batch.path)
    if(!cli)ij.IJ.showMessage('Finished','Images and exact measurements saved in:\n'+batch.path+'\n\nFor wide CSV tables, run scripts/summarize.py (see README).')
    return batch
}
if(cli){try{perform();System.exit(0)}catch(Throwable e){e.printStackTrace();System.exit(1)}}else{perform()}
