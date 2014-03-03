import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

### conditions
from Configuration.AlCa.GlobalTag import *
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# same as input file
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup')
# or get it automatically
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

### standard includes (can't live without)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryAll_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### define tracking sequences to run
# need to recreate hits from clusters (which are stored in RECO)
process.clustToHits = cms.Sequence(
        process.siPixelRecHits*process.siStripMatchedRecHits
        )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/work/r/rovere/TrackingPOG/Studies/CMSSW_7_0_0_pre12/src/reco_trk.root'
#        'file:../../tracks_and_vertices.root'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('trackParameters.root')
                                   )

process.trkproperties = cms.EDAnalyzer('TrackProperties',
                                       TTRHBuilder = cms.string("WithAngleAndTemplate"))


process.p = cms.Path(process.clustToHits
                     * process.trkproperties)
