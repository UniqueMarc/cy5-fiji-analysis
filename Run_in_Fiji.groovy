#@ File(label="This release folder", style="directory") releaseFolder
#@ File(label="Input parent containing group folders", style="directory") inputFolder
#@ File(label="Output parent folder", style="directory") outputFolder
#@ File(label="Configuration JSON (parameters/batch_20260930.json)") configFile
#@ Boolean(label="Reproduce final 20260930 selection using data_manifest.json", value=true) useFinalManifest

// Open in Fiji Script Editor, set language Groovy, and Run.
new GroovyShell(this.class.classLoader,new Binding([
    releaseFolder:releaseFolder,inputFolder:inputFolder,outputFolder:outputFolder,
    configFile:configFile,useFinalManifest:useFinalManifest
])).evaluate(new File(releaseFolder,'scripts/run.groovy'))
