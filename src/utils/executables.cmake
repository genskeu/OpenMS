### the directory name
set(directory source/APPLICATIONS/UTILS)

### list all filenames of the directory here
set(UTILS_executables
AccurateMassSearch
ClusterMassTraces
ClusterMassTracesByPrecursor
CVInspector
DatabaseFilter
DecoyDatabase
DeMeanderize
Digestor
DigestorMotif
ERPairFinder
FeatureFinderSuperHirn
FFEval
FuzzyDiff
IDDecoyProbability
IDExtractor
IDMassAccuracy
IDScoreSwitcher
IDSplitter
LabeledEval
LowMemPeakPickerHiRes
LowMemPeakPickerHiResRandomAccess
MassCalculator
MetaboliteAdductDecharger
MetaboliteSpectralMatcher
MetaProSIP
MRMPairFinder
MSFraggerAdapter
MSSimulator
MultiplexResolver
MzMLSplitter
NovorAdapter
OpenMSInfo
OpenPepXL
OpenPepXLLF
PeakPickerIterative
PSMFeatureExtractor
QCCalculator
QCEmbedder
QCExporter
QCExtractor
QCImporter
QCMerger
QCShrinker
RNADigestor
RNPxl
RNPxlXICFilter
RNPxlSearch
RTEvaluation
SemanticValidator
SequenceCoverageCalculator
SimpleSearchEngine
SiriusAdapter
SpecLibCreator
SpectraSTSearchAdapter
SvmTheoreticalSpectrumGeneratorTrainer
TICCalculator
TransformationEvaluation
XFDR
XMLValidator
)

if(NOT DISABLE_OPENSWATH)
  set(UTILS_executables
    ${UTILS_executables}
    TargetedFileConverter
    OpenSwathDIAPreScoring
    OpenSwathMzMLFileCacher
    OpenSwathWorkflow
    OpenSwathFileSplitter
    OpenSwathRewriteToFeatureXML
    MRMTransitionGroupPicker
  )
endif(NOT DISABLE_OPENSWATH)


## all targets requiring OpenMS_GUI
set(UTILS_executables_with_GUIlib
ImageCreator
INIUpdater
)

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${UTILS_executables} ${UTILS_executables_with_GUIlib})
	list(APPEND sources_VS "${i}.cpp")
endforeach(i)

source_group("" FILES ${sources_VS})
