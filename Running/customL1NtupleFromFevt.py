import FWCore.ParameterSet.Config as cms

from L1TriggerDPG.L1Menu.customL1Ntuple_cfg import *

process.p.remove(process.l1RecoTreeProducer)
#process.p.remove(process.l1MuonRecoTreeProducer)

process.l1MuonRecoTreeProducer.runOnPostLS1 = options.runOnPostLS1
process.l1MuonRecoTreeProducer.triggerMatching = True

# edit here
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles.extend( ['file:///data2/battilan/L1Trigger/62X_RAW_RECO.root'] )

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'ERROR'

