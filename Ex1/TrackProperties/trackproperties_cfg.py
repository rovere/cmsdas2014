import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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

process.trkproperties = cms.EDAnalyzer('TrackProperties'
)


process.p = cms.Path(process.trkproperties)
