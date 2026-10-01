
/* Portable expansion of the executed 20260930 all-visible pipeline.
 * Called by run.groovy. Stage-1 polygon measurements are provisional only.
 * Final measurements come from exact_pixels.groovy; see README.md.
 */
import ij.IJ
import ij.ImagePlus
import ij.ImageStack
import ij.Prefs
import ij.CompositeImage
import ij.gui.Overlay
import ij.gui.Roi
import ij.gui.TextRoi
import ij.io.FileSaver
import ij.measure.Measurements
import ij.measure.ResultsTable
import ij.plugin.RoiEnlarger
import ij.plugin.RoiScaler
import ij.plugin.ZProjector
import ij.plugin.filter.MaximumFinder
import ij.plugin.filter.ParticleAnalyzer
import ij.plugin.frame.RoiManager
import ij.process.ColorProcessor
import ij.process.ImageProcessor
import ij.process.LUT
import loci.plugins.BF

import javax.swing.JFileChooser
import javax.swing.SwingUtilities

import java.awt.Color
import java.awt.Font
import java.awt.GraphicsEnvironment
import java.text.SimpleDateFormat
import java.util.Arrays
import java.util.Date
import java.util.LinkedHashMap
import java.util.HashSet
import java.util.Locale
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream


// Input files are validated and assigned stable IDs by run.groovy.
File[] inputFolders = selectedInputFolders as File[]

// Remove duplicate selections while preserving a reproducible folder order.
Map<String, File> uniqueInputFolders = new LinkedHashMap<String, File>()
for (File folder : inputFolders) {
    if (folder != null && folder.isDirectory()) {
        uniqueInputFolders.put(folder.getCanonicalPath().toLowerCase(Locale.ROOT), folder)
    }
}
inputFolders = uniqueInputFolders.values().toArray(new File[0])
Arrays.sort(inputFolders, { File a, File b ->
    a.getCanonicalPath().compareToIgnoreCase(b.getCanonicalPath())
} as Comparator<File>)
if (inputFolders.length == 0) {
    IJ.error("Input folder error", "None of the selected items is a readable folder.")
    return
}


// ---------- Validate batch settings ----------

if (outputFolder == null) {
    IJ.error("Output folder error", "Choose a parent output folder.")
    return
}
if (hoechstChannel == cy5Channel) {
    IJ.error("Channel error", "Hoechst and Cy5 cannot use the same channel number.")
    return
}
if (maxNucleusAreaUm2 <= minNucleusAreaUm2) {
    IJ.error("Nucleus area error", "Maximum nucleus area must exceed minimum nucleus area.")
    return
}
if (maxPunctumAreaUm2 <= minPunctumAreaUm2) {
    IJ.error("Punctum area error", "Maximum punctum area must exceed minimum punctum area.")
    return
}
if (fixedCy5DisplayMaxRaw <= 0 && fixedCy5DisplayMinRaw > 0) {
    IJ.error("Cy5 display range error",
            "Enter both fixed Cy5 limits, or leave both at 0 for whole-batch automatic calculation.")
    return
}
if (fixedCy5DisplayMaxRaw > 0 &&
        fixedCy5DisplayMaxRaw <= fixedCy5DisplayMinRaw) {
    IJ.error("Cy5 display range error",
            "Fixed Cy5 display maximum must be greater than the minimum.")
    return
}

// ---------- Collect ND2 files from every selected folder ----------

List<Map> nd2Entries = suppliedInputEntries
File[] nd2Files = nd2Entries.collect { (File) it.file }.toArray(new File[0])

String runStamp = new SimpleDateFormat("yyyyMMdd_HHmmss_SSS").format(new Date())
File batchRoot = new File(outputFolder, "Cy5_SingleCell_Fixed640_MultiFolder_Batch_" + runStamp)
if(batchRoot.exists()) throw new IllegalStateException('Refuse to overwrite output '+batchRoot)
binding.setVariable('completedBatchRoot',batchRoot)
Map<String, File> outputDirs = new LinkedHashMap<String, File>()
outputDirs.put("hoechst", new File(batchRoot, "405_Hoechst_ZMAX"))
outputDirs.put("cy5", new File(batchRoot, "640_Cy5_ZMAX"))
outputDirs.put("merged", new File(batchRoot, "Merged_QC"))
outputDirs.put("nucleusMask", new File(batchRoot, "Nucleus_Masks"))
outputDirs.put("punctaMask", new File(batchRoot, "Cy5_Puncta_Masks"))
outputDirs.put("rois", new File(batchRoot, "ROI_Files"))
outputDirs.put("perImage", new File(batchRoot, "Legacy_NOT_for_final_quantitation/Per_Image_Tables"))
outputDirs.put("combined", new File(batchRoot, "Legacy_NOT_for_final_quantitation/Summary"))
outputDirs.put("single405", new File(batchRoot, "Single_Cell_Crops/405_Hoechst"))
outputDirs.put("single640", new File(batchRoot, "Single_Cell_Crops/640_Cy5"))
outputDirs.put("singleMerged", new File(batchRoot, "Single_Cell_Crops/Merged"))
outputDirs.put("singleMarked", new File(batchRoot, "Single_Cell_Crops/Cy5_Marked"))
for (File dir : outputDirs.values()) {
    if (!dir.mkdirs() && !dir.isDirectory()) {
        IJ.error("Output error", "Could not create:\n" + dir.getAbsolutePath())
        return
    }
}


// ---------- Combined result containers ----------

List<Map> allCellRows = []
List<Map> allPunctaRows = []
List<Map> fullFieldRegions=[]
List<Map> fullFieldImages=[]
List<Map> imageSummaryRows = []
List<Map> batchLogRows = []
List<Map> allSingleCellRows = []


// ---------- Small reusable helpers ----------

def safeFileBase = { String name ->
    name.replaceFirst(/(?i)\.nd2$/, "").replaceAll(/[^A-Za-z0-9._-]+/, "_")
}

def closeQuietly = { ImagePlus image ->
    if (image != null) {
        try {
            image.changes = false
            image.close()
        } catch (Throwable ignored) {
        }
    }
}

def histogramQuantile = { int[] histogram, double quantile ->
    long histogramTotal = 0
    for (int count : histogram) histogramTotal += count
    long target = (long) Math.ceil(histogramTotal * quantile)
    long cumulative = 0
    for (int value = 0; value < histogram.length; value++) {
        cumulative += histogram[value]
        if (cumulative >= target) return value
    }
    return histogram.length - 1
}

def longHistogramQuantile = { long[] histogram, double quantile ->
    long histogramTotal = 0
    for (long count : histogram) histogramTotal += count
    if (histogramTotal <= 0) return 0
    long target = (long) Math.ceil(histogramTotal * quantile)
    long cumulative = 0
    for (int value = 0; value < histogram.length; value++) {
        cumulative += histogram[value]
        if (cumulative >= target) return value
    }
    return histogram.length - 1
}

def meanOrNaN = { List<Map> rows, String key ->
    if (rows == null || rows.isEmpty()) return Double.NaN
    double sum = 0.0
    int count = 0
    for (Map row : rows) {
        Object value = row.get(key)
        if (value instanceof Number) {
            double number = ((Number) value).doubleValue()
            if (!Double.isNaN(number) && !Double.isInfinite(number)) {
                sum += number
                count++
            }
        }
    }
    return count > 0 ? sum / count : Double.NaN
}

def makeDisplayBytes = { ImagePlus sourceImage, Double fixedLow, Double fixedHigh ->
    def sourceProcessor = sourceImage.getProcessor()
    int pixelCount = sourceImage.getWidth() * sourceImage.getHeight()
    float[] sortedValues = new float[pixelCount]
    for (int pixelIndex = 0; pixelIndex < pixelCount; pixelIndex++) {
        sortedValues[pixelIndex] = sourceProcessor.getf(pixelIndex)
    }
    Arrays.sort(sortedValues)
    int lowIndex = (int) Math.floor((pixelCount - 1) * 0.00175)
    int highIndex = (int) Math.ceil((pixelCount - 1) * 0.99825)
    boolean useFixedRange = fixedLow != null && fixedHigh != null && fixedHigh > fixedLow
    double displayLow = useFixedRange ? fixedLow : sortedValues[lowIndex]
    double displayHigh = useFixedRange ? fixedHigh : sortedValues[highIndex]
    if (displayHigh <= displayLow) {
        def stats = sourceImage.getStatistics()
        displayLow = stats.min
        displayHigh = stats.max
    }
    if (displayHigh <= displayLow) displayHigh = displayLow + 1.0

    byte[] displayPixels = new byte[pixelCount]
    double scale = 255.0 / (displayHigh - displayLow)
    for (int pixelIndex = 0; pixelIndex < pixelCount; pixelIndex++) {
        double normalized = (sourceProcessor.getf(pixelIndex) - displayLow) * scale
        displayPixels[pixelIndex] = (byte) ((int) Math.round(
                Math.max(0.0, Math.min(255.0, normalized))))
    }
    return [pixels: displayPixels, low: displayLow, high: displayHigh]
}

def cropAndResizeRgb = { ImagePlus sourceImage, int cropX, int cropY,
                         int cropWidth, int cropHeight, int outputPixels,
                         String title ->
    ImageProcessor sourceProcessor = sourceImage.getProcessor()
    sourceProcessor.setRoi(cropX, cropY, cropWidth, cropHeight)
    ImageProcessor croppedProcessor = sourceProcessor.crop()
    sourceProcessor.resetRoi()
    if (croppedProcessor.getWidth() != outputPixels ||
            croppedProcessor.getHeight() != outputPixels) {
        croppedProcessor.setInterpolationMethod(ImageProcessor.BICUBIC)
        croppedProcessor = croppedProcessor.resize(outputPixels, outputPixels, true)
    }
    return new ImagePlus(title, croppedProcessor)
}


// ---------- Determine ONE Cy5 display range for the entire batch ----------

double batchCy5DisplayLow
double batchCy5DisplayHigh
String batchCy5DisplayRangeSource

if (fixedCy5DisplayMaxRaw > 0) {
    batchCy5DisplayLow = fixedCy5DisplayMinRaw
    batchCy5DisplayHigh = fixedCy5DisplayMaxRaw
    batchCy5DisplayRangeSource = "Fixed raw display limits from the selected configuration (display only)"
} else {
    long[] combinedCy5Histogram = new long[65536]
    int prescannedImages = 0
    for (int scanIndex = 0; scanIndex < nd2Files.length; scanIndex++) {
        File scanFile = nd2Files[scanIndex]
        ImagePlus[] scanSeries = null
        ImagePlus scanImp = null
        ImagePlus scanStack = null
        ImagePlus scanProjection = null
        try {
            IJ.showStatus("Measuring common 640 range " + (scanIndex + 1) +
                    "/" + nd2Files.length + ": " + scanFile.getName())
            scanSeries = BF.openImagePlus(scanFile.getAbsolutePath())
            if (scanSeries == null || scanSeries.length == 0) {
                IJ.log("640 pre-scan skipped (could not open): " + scanFile.getName())
                continue
            }
            scanImp = scanSeries[0]
            if (cy5Channel < 1 || cy5Channel > scanImp.getNChannels()) {
                IJ.log("640 pre-scan skipped (channel missing): " + scanFile.getName())
                continue
            }
            if (scanImp.getNSlices() > 1 && scanImp.getNFrames() > 1) {
                IJ.log("640 pre-scan skipped (ambiguous Z/T): " + scanFile.getName())
                continue
            }
            boolean scanFramesAsZ = scanImp.getNSlices() == 1 && scanImp.getNFrames() > 1
            int scanSectionCount = scanFramesAsZ ?
                    scanImp.getNFrames() : scanImp.getNSlices()
            ImageStack stack = new ImageStack(scanImp.getWidth(), scanImp.getHeight())
            for (int section = 1; section <= scanSectionCount; section++) {
                int sourceZ = scanFramesAsZ ? 1 : section
                int sourceT = scanFramesAsZ ? section : 1
                int sourceIndex = scanImp.getStackIndex(cy5Channel, sourceZ, sourceT)
                stack.addSlice(scanImp.getStack().getProcessor(sourceIndex).duplicate())
            }
            scanStack = new ImagePlus("Cy5_batch_range_scan", stack)
            ZProjector scanProjector = new ZProjector(scanStack)
            scanProjector.setMethod(ZProjector.MAX_METHOD)
            scanProjector.setStartSlice(1)
            scanProjector.setStopSlice(scanStack.getStackSize())
            scanProjector.doProjection()
            scanProjection = scanProjector.getProjection()
            int[] histogram = scanProjection.getProcessor().getHistogram()
            int usableBins = Math.min(histogram.length, combinedCy5Histogram.length)
            for (int bin = 0; bin < usableBins; bin++) {
                combinedCy5Histogram[bin] += histogram[bin]
            }
            prescannedImages++
        } catch (Throwable scanError) {
            IJ.log("640 pre-scan skipped: " + scanFile.getName() + " | " +
                    (scanError.getMessage() == null ? scanError.toString() :
                            scanError.getMessage()))
        } finally {
            closeQuietly(scanProjection)
            closeQuietly(scanStack)
            if (scanSeries != null) {
                for (ImagePlus opened : scanSeries) closeQuietly(opened)
            } else {
                closeQuietly(scanImp)
            }
        }
    }
    if (prescannedImages == 0) {
        IJ.error("Cy5 display pre-scan failed",
                "No ND2 image could be used to calculate a common 640 display range.")
        return
    }
    batchCy5DisplayLow = (double) longHistogramQuantile(combinedCy5Histogram, 0.00175)
    batchCy5DisplayHigh = (double) longHistogramQuantile(combinedCy5Histogram, 0.99825)
    if (batchCy5DisplayHigh <= batchCy5DisplayLow) {
        batchCy5DisplayHigh = batchCy5DisplayLow + 1.0
    }
    batchCy5DisplayRangeSource = "Whole-batch pooled Cy5 ZMAX percentiles (0.175%-99.825%)"
}

File fixedRangeFile = new File(outputDirs.combined, "Fixed_640_Display_Range.txt")
fixedRangeFile.text = """Cy5_display_range_source=${batchCy5DisplayRangeSource}
Fixed_Cy5_display_min_raw=${batchCy5DisplayLow}
Fixed_Cy5_display_max_raw=${batchCy5DisplayHigh}
Rule=The same raw minimum and maximum are used for every 640, merged, and single-cell figure in this batch.
Reuse=Enter these two values in the script dialog to reproduce the exact display mapping in another batch.
"""
IJ.log("Fixed batch-wide Cy5 display range: " + batchCy5DisplayLow +
        " - " + batchCy5DisplayHigh + " raw")


// ---------- Process one ND2 file ----------

def processOneImage = { Map inputEntry, int imageIndex ->
    File nd2File = (File) inputEntry.file
    String sourceFolder = (String) inputEntry.rootLabel
    outputDirs.each { key, currentDir ->
 if(key!='combined'){File typeRoot=currentDir.name==~/[0-9]{2}_.*/?currentDir.parentFile:currentDir;File groupDir=new File(typeRoot,sourceFolder);groupDir.mkdirs();outputDirs[key]=groupDir}
 }
 String sourceRelativePath = ((String) inputEntry.relativePath).replace('\\', '/')
    ImagePlus imp = null
    ImagePlus[] openedSeries = null
    ImagePlus hoechstStack = null
    ImagePlus cy5Stack = null
    ImagePlus hoechstRawProjection = null
    ImagePlus cy5Projection = null
    ImagePlus nucleusWork = null
    ImagePlus nucleusMask = null
    ImagePlus cy5Mask = null
    ImagePlus hoechstOnlyImage = null
    ImagePlus cy5OnlyImage = null
    ImagePlus qcComposite = null
    ImagePlus flattenedQc = null
    RoiManager roiManager = null

    String imageName = sourceFolder + "/" + sourceRelativePath
    String fileKey = String.format("%04d_%s", imageIndex, safeFileBase(imageName))

    try {
        IJ.log("\n--- [" + imageIndex + "/" + nd2Files.length + "] " + imageName + " ---")
        IJ.showStatus("Processing " + imageIndex + "/" + nd2Files.length + ": " + imageName)

        def meta=loci.formats.MetadataTools.createOMEXMLMetadata();def reader=new loci.formats.ImageReader();reader.setMetadataStore(meta)
        try{reader.setId(nd2File.path);assert reader.seriesCount==1 && reader.sizeT==1 : 'Unsupported series/time dimensions';assert reader.sizeC>=Math.max(hoechstChannel,cy5Channel) && meta.getChannelName(0,hoechstChannel-1)=='405' && meta.getChannelName(0,cy5Channel-1)=='640' : 'Channel mapping differs';assert reader.pixelType==loci.formats.FormatTools.UINT16 : 'Expected UINT16'}finally{reader.close()}
        openedSeries = BF.openImagePlus(nd2File.getAbsolutePath())
        if (openedSeries == null || openedSeries.length == 0) {
            throw new RuntimeException("Bio-Formats could not open this ND2 file.")
        }
        imp = openedSeries[0]
        if (openedSeries.length > 1) {
            IJ.log("Warning: " + imageName + " contains " + openedSeries.length +
                    " series; only the first series was analyzed.")
        }
        if (imp.getNSlices() > 1 && imp.getNFrames() > 1) {
            throw new RuntimeException("Both Z=" + imp.getNSlices() + " and T=" +
                    imp.getNFrames() + " are greater than 1; optical-section axis is ambiguous.")
        }
        if (hoechstChannel < 1 || hoechstChannel > imp.getNChannels() ||
                cy5Channel < 1 || cy5Channel > imp.getNChannels()) {
            throw new RuntimeException("Image has " + imp.getNChannels() +
                    " channels; requested Hoechst=" + hoechstChannel +
                    " and Cy5=" + cy5Channel + ".")
        }

        def cal = imp.getCalibration()
        if (cal.pixelWidth <= 0 || cal.pixelHeight <= 0 || cal.pixelDepth <= 0) {
            throw new RuntimeException("ND2 spatial calibration is missing or invalid.")
        }
        if(!(cal.getUnit() in ['micron','microns','um','µm','μm'])) throw new IllegalStateException('Expected micrometre XY calibration, got '+cal.getUnit())
        double xyPixelUm = (cal.pixelWidth + cal.pixelHeight) / 2.0
        double pixelAreaUm2 = cal.pixelWidth * cal.pixelHeight
        int width = imp.getWidth()
        int height = imp.getHeight()

        // Force the 640/Cy5 display to magenta for consistent figures.
        LUT cy5DisplayLut = LUT.createLutFromColor(new Color(255, 0, 255))

        boolean sectionsStoredAsFrames = imp.getNSlices() == 1 && imp.getNFrames() > 1
        int opticalSectionCount = sectionsStoredAsFrames ? imp.getNFrames() : imp.getNSlices()
        String opticalSectionSource = sectionsStoredAsFrames ? "T interpreted as Z" : "Z"

        def extractChannelSectionStack = { int channel, String title ->
            ImageStack stack = new ImageStack(width, height)
            for (int section = 1; section <= opticalSectionCount; section++) {
                int sourceZ = sectionsStoredAsFrames ? 1 : section
                int sourceT = sectionsStoredAsFrames ? section : 1
                int stackIndex = imp.getStackIndex(channel, sourceZ, sourceT)
                stack.addSlice(String.format("section_%03d", section),
                        imp.getStack().getProcessor(stackIndex).duplicate())
            }
            ImagePlus extracted = new ImagePlus(title, stack)
            extracted.setCalibration(cal.copy())
            return extracted
        }

        def maximumProjectAllSections = { ImagePlus sectionStack, String title ->
            ZProjector projector = new ZProjector(sectionStack)
            projector.setMethod(ZProjector.MAX_METHOD)
            projector.setStartSlice(1)
            projector.setStopSlice(sectionStack.getStackSize())
            projector.doProjection()
            ImagePlus projection = projector.getProjection()
            projection.setTitle(title)
            projection.setCalibration(cal.copy())
            return projection
        }

        // STEP 1: project every optical section before segmentation/counting.
        hoechstStack = extractChannelSectionStack(hoechstChannel, "Hoechst_all_sections")
        cy5Stack = extractChannelSectionStack(cy5Channel, "Cy5_all_sections")
        hoechstRawProjection = maximumProjectAllSections(hoechstStack,
                "Hoechst_ZMAX_all_" + opticalSectionCount + "_sections")
        cy5Projection = maximumProjectAllSections(cy5Stack,
                "Cy5_ZMAX_all_" + opticalSectionCount + "_sections_for_counting")

        new FileSaver(hoechstRawProjection).saveAsTiff(new File(outputDirs.hoechst,
                fileKey + "_405_Hoechst_ZMAX_RAW16.tif").getAbsolutePath())
        new FileSaver(cy5Projection).saveAsTiff(new File(outputDirs.cy5,
                fileKey + "_640_Cy5_ZMAX_RAW16.tif").getAbsolutePath())

        // STEP 2: segment nuclei on the Hoechst ZMAX.
        nucleusWork = hoechstRawProjection.duplicate()
        nucleusWork.setTitle("Hoechst_nucleus_work")
        double backgroundPx = nucleusBackgroundUm / xyPixelUm
        double sigmaPx = nucleusSigmaUm / xyPixelUm
        if (backgroundPx > 0) {
            IJ.run(nucleusWork, "Subtract Background...", "rolling=" + backgroundPx + " sliding")
        }
        if (sigmaPx > 0) IJ.run(nucleusWork, "Gaussian Blur...", "sigma=" + sigmaPx)

        nucleusMask = nucleusWork.duplicate()
        nucleusMask.setTitle("Hoechst_nucleus_mask")
        def nucleusInputStats = nucleusMask.getStatistics()
        IJ.setAutoThreshold(nucleusMask, nucleusThresholdMethod + " dark")
        double autoRangeMin = nucleusMask.getProcessor().getMinThreshold()
        double autoRangeMax = nucleusMask.getProcessor().getMaxThreshold()
        double autoNucleusThreshold = autoRangeMin <= nucleusInputStats.min + 0.5 ?
                autoRangeMax : autoRangeMin
        double appliedNucleusThreshold = nucleusInputStats.min +
                (autoNucleusThreshold - nucleusInputStats.min) * nucleusThresholdScale
        nucleusMask.getProcessor().setThreshold(appliedNucleusThreshold,
                nucleusInputStats.max, ImageProcessor.NO_LUT_UPDATE)
        Prefs.blackBackground = true
        IJ.run(nucleusMask, "Convert to Mask", "")
        IJ.run(nucleusMask, "Close", "")
        IJ.run(nucleusMask, "Fill Holes", "")
        // Crucial change: do not split one internally uneven nucleus by default.
        if (separateTouchingNuclei) IJ.run(nucleusMask, "Watershed", "")
        nucleusMask.getProcessor().setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
        String nucleusMaskPath = new File(outputDirs.nucleusMask,
                fileKey + "_Hoechst_nucleus_binary_mask.tif").getAbsolutePath()
        new FileSaver(nucleusMask).saveAsTiff(nucleusMaskPath)

        roiManager = RoiManager.getInstance2()
        if (roiManager == null) roiManager = new RoiManager(true)
        roiManager.reset()
        ParticleAnalyzer.setRoiManager(roiManager)
        double minNucleusPixels = minNucleusAreaUm2 / pixelAreaUm2
        double maxNucleusPixels = maxNucleusAreaUm2 / pixelAreaUm2
        ResultsTable nucleusMeasurements = new ResultsTable()
        ParticleAnalyzer nucleusAnalyzer = new ParticleAnalyzer(
                ParticleAnalyzer.ADD_TO_MANAGER |
                        ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES |
                        ParticleAnalyzer.SHOW_NONE,
                Measurements.AREA | Measurements.CENTROID,
                nucleusMeasurements, minNucleusPixels, maxNucleusPixels, 0.10, 1.00)
        nucleusAnalyzer.analyze(nucleusMask)

        int nCells = roiManager.getCount()
        if(nCells==0) IJ.log("No accepted nuclei; export full-field regions, no cell crops.")

        List<Roi> nucleusRois = []
        List<Roi> analysisRois = []
        double[] nucleusCxPx = new double[nCells]
        double[] nucleusCyPx = new double[nCells]
        double[] nucleusArea = new double[nCells]
        double nucleusPaddingPx = nucleusPaddingUm / xyPixelUm
        double expansionPx = expansionUm / xyPixelUm

        for (int i = 0; i < nCells; i++) {
            Roi segmented = (Roi) roiManager.getRoi(i).clone()
            Roi nucleus = nucleusPaddingPx > 0 ?
                    RoiEnlarger.enlarge(segmented, nucleusPaddingPx) : segmented
            if (nucleus == null) nucleus = segmented
            nucleus.setName(String.format("Cell_%03d_nucleus", i + 1))
            nucleusRois.add(nucleus)
            double[] centroid = nucleus.getContourCentroid()
            nucleusCxPx[i] = centroid[0]
            nucleusCyPx[i] = centroid[1]
            nucleusArea[i] = nucleus.getStatistics().area * pixelAreaUm2

            Roi expanded = expansionPx > 0 ? RoiEnlarger.enlarge(nucleus, expansionPx) :
                    (Roi) nucleus.clone()
            if (expanded == null) expanded = (Roi) nucleus.clone()
            expanded.setName(String.format("Cell_%03d_perinuclear_region", i + 1))
            analysisRois.add(expanded)
        }

        // Build mutually exclusive cell territories. Nearest nucleus wins overlap.
        int totalPixels = width * height
        int[] territoryOwner = new int[totalPixels]
        float[] territoryBestDistance2 = new float[totalPixels]
        Arrays.fill(territoryOwner, -1)
        Arrays.fill(territoryBestDistance2, Float.POSITIVE_INFINITY)
        long[] territoryPixelCount = new long[nCells]

        for (int i = 0; i < nCells; i++) {
            Roi region = analysisRois[i]
            Roi nucleus = nucleusRois[i]
            def bounds = region.getBounds()
            int xStart = Math.max(0, bounds.x)
            int yStart = Math.max(0, bounds.y)
            int xEnd = Math.min(width, bounds.x + bounds.width)
            int yEnd = Math.min(height, bounds.y + bounds.height)
            for (int y = yStart; y < yEnd; y++) {
                for (int x = xStart; x < xEnd; x++) {
                    if (!region.contains(x, y)) continue
                    if (excludeNuclearInterior && nucleus.contains(x, y)) continue
                    int pixelIndex = y * width + x
                    double dx = x - nucleusCxPx[i]
                    double dy = y - nucleusCyPx[i]
                    float distance2 = (float) (dx * dx + dy * dy)
                    if (distance2 < territoryBestDistance2[pixelIndex]) {
                        int previousOwner = territoryOwner[pixelIndex]
                        if (previousOwner >= 0) territoryPixelCount[previousOwner]--
                        territoryOwner[pixelIndex] = i
                        territoryBestDistance2[pixelIndex] = distance2
                        territoryPixelCount[i]++
                    }
                }
            }
        }

        // STEP 3: detect puncta on the same Cy5 ZMAX that is exported.
        int[] cy5Histogram = cy5Projection.getProcessor().getHistogram()
        int cy5Median = (int) histogramQuantile(cy5Histogram, 0.5)
        int[] cy5MadHistogram = new int[cy5Histogram.length]
        for (int intensity = 0; intensity < cy5Histogram.length; intensity++) {
            int count = cy5Histogram[intensity]
            if (count > 0) cy5MadHistogram[Math.abs(intensity - cy5Median)] += count
        }
        int cy5Mad = (int) histogramQuantile(cy5MadHistogram, 0.5)
        double cy5RobustSigma = 1.4826 * Math.max(1, cy5Mad)
        double cy5AutoThreshold = cy5Median + cy5ThresholdK * cy5RobustSigma
        double cy5AppliedThreshold = cy5ManualThreshold > 0 ?
                cy5ManualThreshold : cy5AutoThreshold
        def cy5InputStats = cy5Projection.getStatistics()
        // Keep threshold above the maximum when no signal is detectable; do not force a positive.

        def detectionProcessor=cy5Projection.processor.convertToFloatProcessor()
        new ij.plugin.filter.GaussianBlur().blurGaussian(detectionProcessor,cy5DetectionSigmaPixels,cy5DetectionSigmaPixels,0.01d)
        byte[] detectionPixels=new byte[width*height]
        for(int j=0;j<detectionPixels.length;j++)if(detectionProcessor.getf(j)>=cy5AppliedThreshold)detectionPixels[j]=(byte)255
        cy5Mask=new ImagePlus('All_Cy5_connected_regions',new ij.process.ByteProcessor(width,height,detectionPixels,null));cy5Mask.calibration=cal.copy()
        Prefs.blackBackground = true
        cy5Mask.getProcessor().setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
        new FileSaver(cy5Mask).saveAsTiff(new File(outputDirs.punctaMask,
                fileKey + "_Cy5_ZMAX_puncta_mask.tif").getAbsolutePath())

        roiManager.reset()
        ParticleAnalyzer.setRoiManager(roiManager)
        double minPunctumPixels = minPunctumAreaUm2 / pixelAreaUm2
        double maxPunctumPixels = maxPunctumAreaUm2 / pixelAreaUm2
        ResultsTable punctaMeasurements = new ResultsTable()
        ParticleAnalyzer punctaAnalyzer = new ParticleAnalyzer(
                ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.SHOW_NONE,
                Measurements.AREA | Measurements.CENTROID,
                punctaMeasurements, minPunctumPixels, maxPunctumPixels, 0.0, 1.0)
        punctaAnalyzer.analyze(cy5Mask)

        List<Roi> detectedPunctaRois = []
        for (int i = 0; i < roiManager.getCount(); i++) {
            detectedPunctaRois.add((Roi) roiManager.getRoi(i).clone())
        }

        List<Map> thisFullFieldRegions=[]
        int regionId=0
        for(Roi r:detectedPunctaRois){
          regionId++;cy5Projection.setRoi(r);def st=cy5Projection.getStatistics(Measurements.MEAN|Measurements.MIN_MAX);double pixels=r.getStatistics().area;def xy=r.getContourCentroid()
          thisFullFieldRegions.add([Group:inputEntry.group,Source_folder:sourceFolder,Image_Index:imageIndex,Image:fileKey,Source_path:nd2File.path,Region_ID:regionId,Unique_ID:String.format('%04d_R%04d',imageIndex,regionId),Area_um2:pixels*pixelAreaUm2,Area_pixels:(int)pixels,Mean_raw_intensity:st.mean,Max_raw_intensity:st.max,Integrated_density_raw:st.mean*pixels,X_um:xy[0]*cal.pixelWidth,Y_um:xy[1]*cal.pixelHeight])
        }
        cy5Projection.deleteRoi()
        int[] cellCounts = new int[nCells]
        List<Map> assignedSpots = []
        int totalDetected = detectedPunctaRois.size()
        int totalAssigned = 0
        for (Roi detectedRoi : detectedPunctaRois) {
            double[] centroid = detectedRoi.getContourCentroid()
            double xPx = centroid[0]
            double yPx = centroid[1]
            int xi = (int) Math.round(xPx)
            int yi = (int) Math.round(yPx)
            int owner = -1
            if (xi >= 0 && xi < width && yi >= 0 && yi < height) {
                owner = territoryOwner[yi * width + xi]
            }
            if (owner >= 0) {
                totalAssigned++
                cellCounts[owner]++
                Roi assignedRoi = (Roi) detectedRoi.clone()
                assignedRoi.setName(String.format("Cell_%03d_punctum_%04d",
                        owner + 1, totalAssigned))
                cy5Projection.setRoi(assignedRoi)
                def punctumStats = cy5Projection.getStatistics(
                        Measurements.MEAN | Measurements.MIN_MAX)
                double areaPixels = assignedRoi.getStatistics().area
                assignedSpots.add([
                        cell: owner,
                        xUm: (xPx - cal.xOrigin) * cal.pixelWidth,
                        yUm: (yPx - cal.yOrigin) * cal.pixelHeight,
                        roi: assignedRoi,
                        areaUm2: areaPixels * pixelAreaUm2,
                        meanIntensity: punctumStats.mean,
                        maxIntensity: punctumStats.max,
                        integratedDensity: punctumStats.mean * areaPixels
                ])
            }
        }
        cy5Projection.deleteRoi()

        // Per-image CSVs and combined flat rows.
        ResultsTable cellTable = new ResultsTable()
        for (int i = 0; i < nCells; i++) {
            double territoryAreaUm2 = territoryPixelCount[i] * pixelAreaUm2
            double density = territoryAreaUm2 > 0 ?
                    cellCounts[i] * 100.0 / territoryAreaUm2 : Double.NaN
            Map cellRow = [
                    Source_folder: sourceFolder,
                    Relative_path: sourceRelativePath,
                    Image: imageName,
                    Image_Index: imageIndex,
                    Detected_nuclei: nCells,
                    Cell_ID: i + 1,
                    Nucleus_area_um2: nucleusArea[i],
                    Analysis_territory_area_um2: territoryAreaUm2,
                    Cy5_puncta_count: cellCounts[i],
                    Cy5_puncta_per_100um2: density,
                    Nucleus_centroid_X_um: (nucleusCxPx[i] - cal.xOrigin) * cal.pixelWidth,
                    Nucleus_centroid_Y_um: (nucleusCyPx[i] - cal.yOrigin) * cal.pixelHeight
            ]
            allCellRows.add(cellRow)
            cellTable.incrementCounter()
            for (entry in cellRow.entrySet()) {
                if (entry.value instanceof Number) {
                    cellTable.addValue(entry.key.toString(), ((Number) entry.value).doubleValue())
                } else {
                    cellTable.addValue(entry.key.toString(), entry.value.toString())
                }
            }
        }
        cellTable.save(new File(outputDirs.perImage,
                fileKey + "_cell_results.csv").getAbsolutePath())

        ResultsTable punctaTable = new ResultsTable()
        int punctumId = 0
        List<Map> thisImagePunctaRows = []
        for (Map record : assignedSpots) {
            punctumId++
            Map punctaRow = [
                    Source_folder: sourceFolder,
                    Relative_path: sourceRelativePath,
                    Image: imageName,
                    Image_Index: imageIndex,
                    Detected_nuclei: nCells,
                    Punctum_ID: punctumId,
                    Cell_ID: ((int) record.cell) + 1,
                    X_um: (double) record.xUm,
                    Y_um: (double) record.yUm,
                    Detection_image: "Cy5 maximum-intensity Z projection",
                    Area_um2: (double) record.areaUm2,
                    Mean_raw_intensity: (double) record.meanIntensity,
                    Max_raw_intensity: (double) record.maxIntensity,
                    Integrated_density_raw: (double) record.integratedDensity
            ]
            allPunctaRows.add(punctaRow)
            thisImagePunctaRows.add(punctaRow)
            punctaTable.incrementCounter()
            for (entry in punctaRow.entrySet()) {
                if (entry.value instanceof Number) {
                    punctaTable.addValue(entry.key.toString(), ((Number) entry.value).doubleValue())
                } else {
                    punctaTable.addValue(entry.key.toString(), entry.value.toString())
                }
            }
        }
        punctaTable.save(new File(outputDirs.perImage,
                fileKey + "_puncta_results.csv").getAbsolutePath())

        // ROI files.
        roiManager.reset()
        for (Roi roi : nucleusRois) roiManager.addRoi(roi)
        roiManager.runCommand("Save", new File(outputDirs.rois,
                fileKey + "_nucleus_ROIs.zip").getAbsolutePath())
        roiManager.reset()
        for (Roi roi : analysisRois) roiManager.addRoi(roi)
        roiManager.runCommand("Save", new File(outputDirs.rois,
                fileKey + "_perinuclear_analysis_ROIs.zip").getAbsolutePath())

        // Export every accepted particle, including those not assigned to a
        // nucleus. The exact-pixel stage requires the same stable scan order.
        def allRoiZip = new ZipOutputStream(new FileOutputStream(new File(outputDirs.rois,
                fileKey + "_ALL_Cy5_ROIs.zip")))
        try {
            detectedPunctaRois.eachWithIndex { Roi r, int j ->
                allRoiZip.putNextEntry(new ZipEntry(String.format('%04d_R%04d.roi',imageIndex,j+1)))
                new ij.io.RoiEncoder(allRoiZip).write(r)
                allRoiZip.closeEntry()
            }
        } finally { allRoiZip.close() }

        // STEP 4: clean 405, clean 640, unmarked merge, and marked QC.
        Map hoechstDisplay = makeDisplayBytes(hoechstRawProjection, null, null)
        Map cy5Display = makeDisplayBytes(cy5Projection,
                batchCy5DisplayLow, batchCy5DisplayHigh)
        byte[] hoechstIntensity = (byte[]) hoechstDisplay.pixels
        byte[] cy5Intensity = (byte[]) cy5Display.pixels
        byte[] zeroPixels = new byte[hoechstIntensity.length]

        ColorProcessor hoechstProcessor = new ColorProcessor(width, height)
        hoechstProcessor.setRGB(zeroPixels, zeroPixels,
                (byte[]) hoechstIntensity.clone())
        hoechstOnlyImage = new ImagePlus(fileKey + "_405_Hoechst_ZMAX", hoechstProcessor)
        new FileSaver(hoechstOnlyImage).saveAsPng(new File(outputDirs.hoechst,
                fileKey + "_405_Hoechst_ZMAX.png").getAbsolutePath())

        byte[] cy5OnlyRed = new byte[cy5Intensity.length]
        byte[] cy5OnlyGreen = new byte[cy5Intensity.length]
        byte[] cy5OnlyBlue = new byte[cy5Intensity.length]
        byte[] mergeRed = new byte[cy5Intensity.length]
        byte[] mergeGreen = new byte[cy5Intensity.length]
        byte[] mergeBlue = new byte[cy5Intensity.length]
        for (int p = 0; p < cy5Intensity.length; p++) {
            int cy5Value = cy5Intensity[p] & 0xff
            int hoechstValue = hoechstIntensity[p] & 0xff
            int r = cy5DisplayLut.getRed(cy5Value)
            int g = cy5DisplayLut.getGreen(cy5Value)
            int b = cy5DisplayLut.getBlue(cy5Value)
            cy5OnlyRed[p] = (byte) r
            cy5OnlyGreen[p] = (byte) g
            cy5OnlyBlue[p] = (byte) b
            mergeRed[p] = (byte) r
            mergeGreen[p] = (byte) g
            mergeBlue[p] = (byte) Math.min(255, hoechstValue + b)
        }
        ColorProcessor cy5OnlyProcessor = new ColorProcessor(width, height)
        cy5OnlyProcessor.setRGB(cy5OnlyRed, cy5OnlyGreen, cy5OnlyBlue)
        cy5OnlyImage = new ImagePlus(fileKey + "_640_Cy5_ZMAX", cy5OnlyProcessor)
        new FileSaver(cy5OnlyImage).saveAsPng(new File(outputDirs.cy5,
                fileKey + "_640_Cy5_ZMAX.png").getAbsolutePath())

        ColorProcessor mergeProcessor = new ColorProcessor(width, height)
        mergeProcessor.setRGB(mergeRed, mergeGreen, mergeBlue)
        qcComposite = new ImagePlus(fileKey + "_ZMAX_merge", mergeProcessor)
        new FileSaver(qcComposite).saveAsPng(new File(outputDirs.merged,
                fileKey + "_ZMAX_merged_unmarked.png").getAbsolutePath())

        // Full-field QC: cyan nuclei and green puncta only. No numbers, total
        // nucleus count, or yellow perinuclear-region outlines.
        Overlay overlay = new Overlay()
        for (int i = 0; i < nCells; i++) {
            Roi nucleusOutline = (Roi) nucleusRois[i].clone()
            nucleusOutline.setStrokeColor(Color.CYAN)
            nucleusOutline.setStrokeWidth(1.0)
            overlay.add(nucleusOutline)
        }
        for (Roi detectedRoi : detectedPunctaRois) {
            Roi marker = (Roi) detectedRoi.clone()
            marker.setStrokeColor(Color.GREEN)
            marker.setStrokeWidth(1.0)
            overlay.add(marker)
        }
        qcComposite.setOverlay(overlay)
        flattenedQc = qcComposite.flatten()
        new FileSaver(flattenedQc).saveAsPng(new File(outputDirs.merged,
                fileKey + "_QC_ZMAX_marked.png").getAbsolutePath())

        // STEP 5: export standardized single-cell fields.
        int cropWidthPx = Math.max(2,
                (int) Math.round(singleCellCropUm / cal.pixelWidth))
        int cropHeightPx = Math.max(2,
                (int) Math.round(singleCellCropUm / cal.pixelHeight))
        double cropScaleX = singleCellOutputPixels / (double) cropWidthPx
        double cropScaleY = singleCellOutputPixels / (double) cropHeightPx
        int exportedSingleCells = 0

        for (int cellIndex = 0; cellIndex < nCells; cellIndex++) {
            double nearestNeighborUm = Double.POSITIVE_INFINITY
            for (int otherIndex = 0; otherIndex < nCells; otherIndex++) {
                if (otherIndex == cellIndex) continue
                double dxUm = (nucleusCxPx[cellIndex] - nucleusCxPx[otherIndex]) *
                        cal.pixelWidth
                double dyUm = (nucleusCyPx[cellIndex] - nucleusCyPx[otherIndex]) *
                        cal.pixelHeight
                double distanceUm = Math.sqrt(dxUm * dxUm + dyUm * dyUm)
                if (distanceUm < nearestNeighborUm) nearestNeighborUm = distanceUm
            }
            if (nearestNeighborUm < minimumNeighborDistanceUm) continue

            int cropX = (int) Math.round(nucleusCxPx[cellIndex] - cropWidthPx / 2.0)
            int cropY = (int) Math.round(nucleusCyPx[cellIndex] - cropHeightPx / 2.0)
            // A standardized crop must be complete; never pad missing image edges.
            if (cropX < 0 || cropY < 0 || cropX + cropWidthPx > width ||
                    cropY + cropHeightPx > height) continue

            String cellKey = fileKey + String.format("_Cell_%03d", cellIndex + 1)
            ImagePlus cell405 = null
            ImagePlus cell640 = null
            ImagePlus cellMerged = null
            ImagePlus cellMarked = null
            try {
                cell405 = cropAndResizeRgb(hoechstOnlyImage, cropX, cropY,
                        cropWidthPx, cropHeightPx, singleCellOutputPixels,
                        cellKey + "_405")
                cell640 = cropAndResizeRgb(cy5OnlyImage, cropX, cropY,
                        cropWidthPx, cropHeightPx, singleCellOutputPixels,
                        cellKey + "_640")
                cellMerged = cropAndResizeRgb(qcComposite, cropX, cropY,
                        cropWidthPx, cropHeightPx, singleCellOutputPixels,
                        cellKey + "_merged")

                new FileSaver(cell405).saveAsPng(new File(outputDirs.single405,
                        cellKey + "_405_Hoechst.png").getAbsolutePath())
                new FileSaver(cell640).saveAsPng(new File(outputDirs.single640,
                        cellKey + "_640_Cy5_magenta.png").getAbsolutePath())
                new FileSaver(cellMerged).saveAsPng(new File(outputDirs.singleMerged,
                        cellKey + "_405_640_merged.png").getAbsolutePath())

                Overlay cellOverlay = new Overlay()
                Roi localNucleus = (Roi) nucleusRois[cellIndex].clone()
                def nucleusBounds = localNucleus.getBounds()
                localNucleus.setLocation(nucleusBounds.x - cropX,
                        nucleusBounds.y - cropY)
                Roi scaledNucleus = RoiScaler.scale(localNucleus,
                        cropScaleX, cropScaleY, false)
                scaledNucleus.setStrokeColor(Color.CYAN)
                scaledNucleus.setStrokeWidth(2.0)
                cellOverlay.add(scaledNucleus)

                for (Roi detectedRoi : detectedPunctaRois) {
                    if(!detectedRoi.bounds.intersects(cropX,cropY,cropWidthPx,cropHeightPx))continue
                    Roi localPunctum = (Roi) detectedRoi.clone()
                    def punctumBounds = localPunctum.getBounds()
                    localPunctum.setLocation(punctumBounds.x - cropX,
                            punctumBounds.y - cropY)
                    Roi scaledPunctum = RoiScaler.scale(localPunctum,
                            cropScaleX, cropScaleY, false)
                    scaledPunctum.setStrokeColor(Color.GREEN)
                    scaledPunctum.setStrokeWidth(2.0)
                    cellOverlay.add(scaledPunctum)
                }
                cellMerged.setOverlay(cellOverlay)
                cellMarked = cellMerged.flatten()
                new FileSaver(cellMarked).saveAsPng(new File(outputDirs.singleMarked,
                        cellKey + "_Cy5_marked_nucleus_outline.png").getAbsolutePath())

                exportedSingleCells++
                allSingleCellRows.add([
                        Source_folder: sourceFolder,
                        Relative_path: sourceRelativePath,
                        Image: imageName,
                        Image_Index: imageIndex,
                        Source_Cell_ID: cellIndex + 1,
                        Crop_ID: cellKey,
                        Center_X_um: (nucleusCxPx[cellIndex] - cal.xOrigin) *
                                cal.pixelWidth,
                        Center_Y_um: (nucleusCyPx[cellIndex] - cal.yOrigin) *
                                cal.pixelHeight,
                        Nearest_neighbor_center_distance_um:
                                Double.isInfinite(nearestNeighborUm) ? "" : nearestNeighborUm,
                        Crop_width_um: singleCellCropUm,
                        Crop_height_um: singleCellCropUm,
                        Output_width_pixels: singleCellOutputPixels,
                        Output_height_pixels: singleCellOutputPixels,
                        Assigned_Cy5_puncta: cellCounts[cellIndex]
                ])
            } finally {
                closeQuietly(cellMarked)
                if (cellMerged != null) cellMerged.setOverlay(null)
                closeQuietly(cellMerged)
                closeQuietly(cell640)
                closeQuietly(cell405)
            }
        }
        IJ.log("Standardized single-cell crops exported: " + exportedSingleCells)

        double meanArea = meanOrNaN(thisImagePunctaRows, "Area_um2")
        double meanIntegrated = meanOrNaN(thisImagePunctaRows, "Integrated_density_raw")
        Map summaryRow = [
                Source_folder: sourceFolder,
                Relative_path: sourceRelativePath,
                Image_Index: imageIndex,
                Image: imageName,
                Status: "OK",
                Error_Message: "",
                Detected_nuclei: nCells,
                Optical_sections_projected: opticalSectionCount,
                Total_Cy5_objects_detected: totalDetected,
                Assigned_puncta: totalAssigned,
                Unassigned_puncta: totalDetected - totalAssigned,
                Single_cell_crops_exported: exportedSingleCells,
                Mean_Area_um2: meanArea,
                Mean_Integrated_density_raw: meanIntegrated
        ]
        imageSummaryRows.add(summaryRow)
        batchLogRows.add([
                Source_folder: sourceFolder,
                Relative_path: sourceRelativePath,
                Image_Index: imageIndex,
                Image: imageName,
                Status: "OK",
                Message: "Completed",
                Detected_nuclei: nCells,
                Assigned_puncta: totalAssigned,
                Single_cell_crops_exported: exportedSingleCells
        ])

        File parameterFile = new File(outputDirs.perImage,
                fileKey + "_parameters.txt")
        parameterFile.text = """Image=${imageName}
Source_path=${nd2File.getAbsolutePath()}
Hoechst_channel=${hoechstChannel}
Cy5_channel=${cy5Channel}
Cy5_display_range_source=${batchCy5DisplayRangeSource}
Fixed_Cy5_display_min_raw=${batchCy5DisplayLow}
Fixed_Cy5_display_max_raw=${batchCy5DisplayHigh}
Original_C=${imp.getNChannels()}
Original_Z=${imp.getNSlices()}
Original_T=${imp.getNFrames()}
Optical_section_source=${opticalSectionSource}
Optical_sections_projected=${opticalSectionCount}
Projection=Maximum intensity across all optical sections before segmentation/counting
Cell_boundary_mode=Fixed-distance perinuclear region; brightfield not used
Overlapping_region_rule=Nearest nucleus owns each pixel; no duplicate punctum assignment
Expansion_um=${expansionUm}
Exclude_nuclear_interior=${excludeNuclearInterior}
Nucleus_threshold_method=${nucleusThresholdMethod}
Nucleus_auto_threshold=${autoNucleusThreshold}
Nucleus_threshold_scale=${nucleusThresholdScale}
Nucleus_applied_threshold=${appliedNucleusThreshold}
Nucleus_area_um2=${minNucleusAreaUm2}-${maxNucleusAreaUm2}
Hoechst_background_radius_um=${nucleusBackgroundUm}
Hoechst_gaussian_sigma_um=${nucleusSigmaUm}
Nucleus_ROI_outward_padding_um=${nucleusPaddingUm}
Nucleus_watershed_enabled=${separateTouchingNuclei}
Cy5_background_median_raw=${cy5Median}
Cy5_background_MAD_raw=${cy5Mad}
Cy5_robust_sigma_raw=${cy5RobustSigma}
Cy5_threshold_K=${cy5ThresholdK}
Cy5_auto_threshold_raw=${cy5AutoThreshold}
Cy5_manual_threshold_raw=${cy5ManualThreshold}
Cy5_applied_threshold_raw=${cy5AppliedThreshold}
Cy5_split_connected_objects=${splitConnectedCy5}
Cy5_peak_separation_tolerance_raw=${cy5PeakTolerance}
Cy5_punctum_area_um2=${minPunctumAreaUm2}-${maxPunctumAreaUm2}
Detected_nuclei=${nCells}
Total_2D_projection_spots_detected=${totalDetected}
Spots_assigned_to_cells=${totalAssigned}
Unassigned_spots=${totalDetected - totalAssigned}
Mean_Area_um2=${meanArea}
Mean_Integrated_density_raw=${meanIntegrated}
Integrated_density_definition=Mean raw Cy5 intensity multiplied by punctum area in pixels
Single_cell_crop_um=${singleCellCropUm}
Single_cell_output_pixels=${singleCellOutputPixels}x${singleCellOutputPixels}
Minimum_neighbor_nucleus_center_distance_um=${minimumNeighborDistanceUm}
Single_cell_crops_exported=${exportedSingleCells}
"""

        fullFieldRegions.addAll(thisFullFieldRegions)
        fullFieldImages.add([Group:inputEntry.group,Source_folder:sourceFolder,Image_Index:imageIndex,Image:fileKey,Source_path:nd2File.path,Regions:totalDetected,Nuclei:nCells,Crops:exportedSingleCells,Z:opticalSectionCount,Pixel_width_um:cal.pixelWidth,Pixel_height_um:cal.pixelHeight,Threshold_raw:cy5AppliedThreshold,Fixed_640_min:batchCy5DisplayLow,Fixed_640_max:batchCy5DisplayHigh])
        new File(batchRoot,'All_field_checkpoint.json').setText(groovy.json.JsonOutput.toJson([regions:fullFieldRegions,images:fullFieldImages]),'UTF-8')
        println('DONE '+imageIndex+'/'+nd2Files.length+' '+sourceFolder+' nuclei='+nCells+' regions='+totalDetected+' crops='+exportedSingleCells)
        IJ.log("Completed: nuclei=" + nCells + ", assigned puncta=" + totalAssigned +
                ", Z sections=" + opticalSectionCount)
        return summaryRow

    } finally {
        if (roiManager != null) {
            try { roiManager.reset() } catch (Throwable ignored) {}
        }
        closeQuietly(flattenedQc)
        closeQuietly(qcComposite)
        closeQuietly(hoechstOnlyImage)
        closeQuietly(cy5OnlyImage)
        closeQuietly(cy5Mask)
        closeQuietly(nucleusMask)
        closeQuietly(nucleusWork)
        closeQuietly(cy5Projection)
        closeQuietly(hoechstRawProjection)
        closeQuietly(cy5Stack)
        closeQuietly(hoechstStack)
        if (openedSeries != null) {
            for (ImagePlus opened : openedSeries) closeQuietly(opened)
        } else {
            closeQuietly(imp)
        }
        System.gc()
    }
}


// ---------- Run every file, continuing after per-image failures ----------

for (int fileIndex = 0; fileIndex < nd2Entries.size(); fileIndex++) {
    Map inputEntry = nd2Entries[fileIndex]
    File nd2File = (File) inputEntry.file
    int imageIndex = ((Number) inputEntry.stableId).intValue()
    try {
        processOneImage(inputEntry, imageIndex)
    } catch (Throwable error) {
        String message = error.getMessage() == null ? error.toString() : error.getMessage()
        String sourceFolder = (String) inputEntry.rootLabel
        outputDirs.each { key, currentDir ->
 if(key!='combined'){File typeRoot=currentDir.name==~/[0-9]{2}_.*/?currentDir.parentFile:currentDir;File groupDir=new File(typeRoot,sourceFolder);groupDir.mkdirs();outputDirs[key]=groupDir}
 }
 String sourceRelativePath = ((String) inputEntry.relativePath).replace('\\', '/')
        String imageName = sourceFolder + "/" + sourceRelativePath
        IJ.log("FAILED: " + imageName + " | " + message)
        imageSummaryRows.add([
                Source_folder: sourceFolder,
                Relative_path: sourceRelativePath,
                Image_Index: imageIndex,
                Image: imageName,
                Status: "FAILED",
                Error_Message: message,
                Detected_nuclei: "",
                Optical_sections_projected: "",
                Total_Cy5_objects_detected: "",
                Assigned_puncta: "",
                Unassigned_puncta: "",
                Single_cell_crops_exported: "",
                Mean_Area_um2: "",
                Mean_Integrated_density_raw: ""
        ])
        batchLogRows.add([
                Source_folder: sourceFolder,
                Relative_path: sourceRelativePath,
                Image_Index: imageIndex,
                Image: imageName,
                Status: "FAILED",
                Message: message,
                Detected_nuclei: "",
                Assigned_puncta: "",
                Single_cell_crops_exported: ""
        ])
    }
}

// One index row per exported single-cell crop set (four PNG files per row).
ResultsTable singleCellIndexTable = new ResultsTable()
for (Map cropRow : allSingleCellRows) {
    singleCellIndexTable.incrementCounter()
    for (entry in cropRow.entrySet()) {
        if (entry.value instanceof Number) {
            double value = ((Number) entry.value).doubleValue()
            if (!Double.isNaN(value) && !Double.isInfinite(value)) {
                singleCellIndexTable.addValue(entry.key.toString(), value)
            }
        } else if (entry.value != null) {
            singleCellIndexTable.addValue(entry.key.toString(), entry.value.toString())
        }
    }
}
singleCellIndexTable.save(new File(new File(batchRoot, "Single_Cell_Crops"),
        "Single_Cell_Index.csv").getAbsolutePath())


// Add an overall row after every per-image summary row.
int successfulImages = 0
int failedImages = 0
int totalNuclei = 0
int totalDetectedObjects = 0
int totalAssignedPuncta = 0
int totalUnassignedPuncta = 0
int totalSingleCellCrops = 0
for (Map row : imageSummaryRows) {
    if (row.Status == "OK") {
        successfulImages++
        totalNuclei += ((Number) row.Detected_nuclei).intValue()
        totalDetectedObjects += ((Number) row.Total_Cy5_objects_detected).intValue()
        totalAssignedPuncta += ((Number) row.Assigned_puncta).intValue()
        totalUnassignedPuncta += ((Number) row.Unassigned_puncta).intValue()
        totalSingleCellCrops += ((Number) row.Single_cell_crops_exported).intValue()
    } else {
        failedImages++
    }
}
imageSummaryRows.add([
        Source_folder: "",
        Relative_path: "",
        Image_Index: "",
        Image: "OVERALL",
        Status: successfulImages + " completed; " + failedImages + " failed",
        Error_Message: "Means are calculated across all assigned puncta from all successful images.",
        Detected_nuclei: totalNuclei,
        Optical_sections_projected: "",
        Total_Cy5_objects_detected: totalDetectedObjects,
        Assigned_puncta: totalAssignedPuncta,
        Unassigned_puncta: totalUnassignedPuncta,
        Single_cell_crops_exported: totalSingleCellCrops,
        Mean_Area_um2: meanOrNaN(allPunctaRows, "Area_um2"),
        Mean_Integrated_density_raw: meanOrNaN(allPunctaRows, "Integrated_density_raw")
])

// Final portable CSV tables are generated separately from exact pixels.
new File(batchRoot,'All_field_results.json').setText(groovy.json.JsonOutput.toJson([regions:fullFieldRegions,images:fullFieldImages,cells:allCellRows,crops:allSingleCellRows,status:batchLogRows,sourceImages:nd2Files.length]),'UTF-8')
println('BATCH_OUTPUT '+batchRoot.path)

File batchReadme = new File(outputDirs.combined, "Batch_Summary.txt")
String inputFolderList = inputFolders.collect { it.getAbsolutePath() }.join(" | ")
batchReadme.text = """Batch_created=${runStamp}
Input_folders=${inputFolderList}
Selected_input_folder_count=${inputFolders.length}
Include_subfolders=${includeSubfolders}
Output_folder=${batchRoot.getAbsolutePath()}
ND2_files_found=${nd2Files.length}
Images_completed=${successfulImages}
Images_failed=${failedImages}
Total_detected_nuclei=${totalNuclei}
Total_assigned_puncta=${totalAssignedPuncta}
Total_single_cell_crops_exported=${totalSingleCellCrops}
Mean_Area_um2_all_assigned_puncta=${meanOrNaN(allPunctaRows, "Area_um2")}
Mean_Integrated_density_raw_all_assigned_puncta=${meanOrNaN(allPunctaRows, "Integrated_density_raw")}
Nucleus_watershed_enabled=${separateTouchingNuclei}
Duplicate_counting_prevention=Overlapping perinuclear pixels assigned to nearest nucleus
Single_cell_crop_um=${singleCellCropUm}
Single_cell_output_pixels=${singleCellOutputPixels}x${singleCellOutputPixels}
Minimum_neighbor_nucleus_center_distance_um=${minimumNeighborDistanceUm}
Cy5_display_range_source=${batchCy5DisplayRangeSource}
Fixed_Cy5_display_min_raw=${batchCy5DisplayLow}
Fixed_Cy5_display_max_raw=${batchCy5DisplayHigh}
"""

IJ.showStatus("Batch analysis finished")
IJ.log("\n=== Batch finished ===")
IJ.log("Completed images: " + successfulImages)
IJ.log("Failed images: " + failedImages)
IJ.log("Single-cell crop sets exported: " + totalSingleCellCrops)
IJ.log("Fixed Cy5 display range: " + batchCy5DisplayLow + " - " +
        batchCy5DisplayHigh + " raw")
IJ.log("Stage 1 complete. Use Exact_pixel_measurements for final quantitation.")
if (false) {
    IJ.showMessage("Batch finished",
            "Input folders: " + inputFolders.length +
            "\nND2 files: " + nd2Files.length +
            "\nCompleted: " + successfulImages +
            "\nFailed: " + failedImages +
            "\nSingle-cell crop sets: " + totalSingleCellCrops +
            "\nFixed Cy5 display: " + batchCy5DisplayLow + " - " +
            batchCy5DisplayHigh + " raw" +
            "\n\nResults saved in:\n" + batchRoot.getAbsolutePath())
}
if(failedImages>0)throw new RuntimeException("Batch includes "+failedImages+" failed images; review status.")

return batchRoot
