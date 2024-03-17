// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief reconstruction of Sigma^0 resonance
/// \author Sergei Solokhin
/// \since 12/03/2024

#include "Framework/AnalysisTask.h" // needed in any case
#include "Framework/runDataProcessing.h" // needed in any case
#include "Framework/ASoA.h" // For creating columns and tables
#include "Framework/AnalysisDataModel.h" // For extending the standard AOD format
#include "Framework/ASoAHelpers.h" // For Filters

#include "Common/DataModel/EventSelection.h" // For event selection
#include "Common/DataModel/TrackSelectionTables.h" // For DCA info
#include "Common/DataModel/PIDResponse.h" // For PID: expected values, resolutions,etc.

#include "PWGLF/DataModel/LFStrangenessTables.h" // for strangeness analysis
#include "CommonConstants/PhysicsConstants.h" // for masses

#include "TLorentzVector.h"
#include <cmath>

using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions; // for filters

namespace o2::aod {
  namespace conversionphoton {
    DECLARE_SOA_COLUMN(Px, px, float);
    DECLARE_SOA_COLUMN(Py, py, float);
    DECLARE_SOA_COLUMN(Pz, pz, float);
    DECLARE_SOA_COLUMN(E, e, float);
  }

  DECLARE_SOA_TABLE(ConversionPhoton, "AOD", "PCMPHOTON", 
  conversionphoton::Px, conversionphoton::Py, conversionphoton::Pz, conversionphoton::E);
}

namespace o2::aod {
  namespace lambdahyperon {
    DECLARE_SOA_COLUMN(Px, px, float);
    DECLARE_SOA_COLUMN(Py, py, float);
    DECLARE_SOA_COLUMN(Pz, pz, float);
    DECLARE_SOA_COLUMN(E, e, float);
  }

  DECLARE_SOA_TABLE(LambdaHyperon, "AOD", "LAMBDAHYPERON", 
  lambdahyperon::Px, lambdahyperon::Py, lambdahyperon::Pz, lambdahyperon::E);
}

struct ProcessCollisions {
  Configurable<float> zVertexCut{"zVertexCut", 10.0f, "Maximum Primary Vertex Z coordinate [cm]"};

  HistogramRegistry histosCollisions{"histosCollisions", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&) {
    const AxisSpec axisVertexZ{100, -20.0f, 20.0f, "Vertex Z Coordinate [cm]"};

    histosCollisions.add("vertexZ", "Primary Vertex Z Coordinate", kTH1F, {axisVertexZ});
    histosCollisions.add("eventCounter", "Event Counter", kTH1D, {{1, 0, 1, " Accepted Events Count"}});
  }

  void process(Collision const& collision) {
    histosCollisions.fill(HIST("vertexZ"), collision.posZ());
    if (abs(collision.posZ()) < zVertexCut && collision.posZ() != 0.0f) histosCollisions.fill(HIST("eventCounter"), 0);
  }
};

struct ProcessGeneratedEvents {
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<float> maxPhotonPt{"maxPhotonPt", 1.0f, "Maximum Photon Transverse Momentum [GeV]"};
  Configurable<float> maxLambdaPt{"maxLambdaPt", 5.0f, "Maximum Lambda Hyperon Transverse Momentum [GeV]"};
  Configurable<float> maxSigma0Pt{"maxSigma0Pt", 5.0f, "Maximum Sigma^0 Hyperon Transverse Momentum [GeV]"};

  HistogramRegistry histosMC{"histosMC", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&) {
    const AxisSpec axisEta{150, -1.5f, 1.5f, "#eta"};
    const AxisSpec axisPhotonPt{nBinsPt, 0.0f, maxPhotonPt, "p_{T, #gamma} [GeV]"};
    const AxisSpec axisLambdaPt{nBinsPt, 0.0f, maxLambdaPt, "p_{T, #Lambda} [GeV]"};
    const AxisSpec axisSigma0Pt{nBinsPt, 0.0f, maxSigma0Pt, "p_{T, #Sigma^{0}} [GeV]"};

    histosMC.add("mcGenEta", "Pseudorapidity Of Generated Tracks (MC)", kTH1F, {axisEta});
    histosMC.add("mcGenPt", "Transverse Momentum Of Generated Tracks (MC)", kTH1F, {{100, 0.0f, 10.0f, "p_{T} [GeV]"}});
    histosMC.add("mcGenPhotonPt", "Transverse Momentum Of Generated Photons (MC)", kTH1F, {axisPhotonPt});    
    histosMC.add("mcGenLambdaPt", "Transverse Momentum Of Generated #Lambda Hyperons (MC)", kTH1F, {axisLambdaPt});
    histosMC.add("mcGenSigma0Pt", "Transverse Momentum Of Generated #Sigma^{0} Hyperons (MC)", kTH1F, {axisSigma0Pt});
    histosMC.add("mcGenPhotonFromSigma0Pt", "Transverse Momentum Of Generated Photons From #Sigma^{0} Hyperons (MC)", kTH1F, {axisPhotonPt});
    histosMC.add("mcGenLambdaFromSigma0Pt", "Transverse Momentum Of Generated #Lambda Hyperons From #Sigma^{0} Hyperons (MC)", kTH1F, {axisLambdaPt});
  }

  void process(McCollision const&, McParticles const& mcparticles) {
    for (auto const& mcparticle: mcparticles) {
      histosMC.fill(HIST("mcGenEta"), mcparticle.eta());
      histosMC.fill(HIST("mcGenPt"), mcparticle.pt());
      if (abs(mcparticle.pdgCode()) == 3122) {
        histosMC.fill(HIST("mcGenLambdaPt"), mcparticle.pt());
        auto const& mother = mcparticle.mothers_first_as<McParticles>();
        if (abs(mother.pdgCode()) == 3212) histosMC.fill(HIST("mcGenLambdaFromSigma0Pt"), mcparticle.pt());
      }

      if (mcparticle.pdgCode() == 22) {
        histosMC.fill(HIST("mcGenPhotonPt"), mcparticle.pt());
        auto const& mother = mcparticle.mothers_first_as<McParticles>();
        if (abs(mother.pdgCode()) == 3212) histosMC.fill(HIST("mcGenPhotonFromSigma0Pt"), mcparticle.pt());
      }

      if (abs(mcparticle.pdgCode()) == 3212) {
        histosMC.fill(HIST("mcGenSigma0Pt"), mcparticle.pt());
      }
    }
  }
};

struct ProcessConversionPhotons {
  Produces<ConversionPhoton> AddConversionPhoton;

  float electronMass = o2::constants::physics::MassElectron;

  Configurable<float> zVertexCut{"zVertexCut", 10.0f, "Maximum Primary Vertex Z coordinate [cm]"};
  Filter zVertexFilter = (nabs(collision::posZ) < zVertexCut);
  Filter zVertexErrorFilter = (collision::posZ != 0.0f);
  Filter eventSelectionFilter = (evsel::sel8 == true);
  using filteredCollision = Filtered<Join<Collisions, EvSels>>::iterator;

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsMass{"nBinsMass", 100, "N bins in invariant mass histo"};
  Configurable<float> etaCut{"etaCut", 1.2f, "Maximum Pseudorapidity"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1.0f, "DCA V0 Daughters [cm]"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06f, "DCA Pos To PV [cm]"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06f, "DCA Neg To PV [cm]"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.995, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.0f, "V0 Radius [cm]"};

  Filter preV0Filter = (nabs(v0data::dcapostopv) > v0setting_dcapostopv && nabs(v0data::dcanegtopv) > v0setting_dcanegtopv && v0data::dcaV0daughters < v0setting_dcav0dau);
  using filteredV0s = Filtered<Join<V0Datas, McV0Labels>>;

  Configurable<float> maxPhotonMass{"maxPhotonMass", 0.1f, "Maximum Electron-Positron Invariant Mass [GeV]"};
  Configurable<float> maxPhotonPt{"maxPhotonPt", 2.5f, "Maximum Photon Transverse Momentum [GeV]"};
  Configurable<float> maxTpcNSigmaEl{"maxTpcNSigmaEl", 3.0f, "Maximum N_{#sigma_{e}} from TPC signal"};
  Configurable<float> maxTrackDCA{"maxTrackDCA", 1.0f, "Maximum Track DCA [cm]"};
  Configurable<float> maxPhotonAlpha{"maxPhotonAlpha", 0.7f, "Maximum Photon Decay Asymmetry"};
  Configurable<float> maxPhotonQt{"maxPhotonQt", 0.05f, "Maximum Photon Decay q_{T} [GeV]"};

  Filter etaFilter = (nabs(track::eta) < etaCut);
  Filter dcaFilter = (nabs(track::dcaXY) < maxTrackDCA);

  using photonDaughterTracks = Join<Tracks, TracksDCA, pidTPCEl, McTrackLabels>;
  using filteredPhotonDaughterTracks = Filtered<photonDaughterTracks>;

  HistogramRegistry histosConversionPhoton{"histosConversionPhoton", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&) {
    const AxisSpec axisEta{100, -etaCut, etaCut, "#eta"};
    const AxisSpec axisPhotonMass{nBinsMass, 0.0f, maxPhotonMass, "M_{e^{+}e^{-}} [GeV]"};
    const AxisSpec axisPhotonPt{nBinsPt, 0.0f, maxPhotonPt, "p_{T} [GeV]"};
    const AxisSpec axisNsigmaE{100, -maxTpcNSigmaEl, maxTpcNSigmaEl, "N_{#sigma_{e}}"};
    const AxisSpec axisDCA{100, -5.0f, 5.0f, "DCA [cm]"};
    const AxisSpec axisV0Radius{100, 0.0f, 50.0f, "V0 Radius [cm]"};
    const AxisSpec axisAlpha{100, -1.0f, 1.0f, "#alpha"};
    const AxisSpec axisQt{100, 0.0f, 0.25f, "q_{T} [GeV]"};

    histosConversionPhoton.add("photonCutEffects", "Cuts Effect On Number Of Kept Photon Candidates", kTH1D, {{10, 0, 10, "Cut"}});

    histosConversionPhoton.add("photonEta", "Pseudorapidity Of Photon Candidates", kTH1F, {axisEta});
    histosConversionPhoton.add("photonMass", "Invariant Mass Of Photon Candidates", kTH1F, {axisPhotonMass});
    histosConversionPhoton.add("photonPt", "Transverse Momentum Of Photon Candidates", kTH1F, {axisPhotonPt});

    histosConversionPhoton.add("posTPCEl", "N_{#sigma, e} For Positively Charged Decay Vertex Daughter Tracks In TPC", kTH1F, {axisNsigmaE});
    histosConversionPhoton.add("negTPCEl", "N_{#sigma, e} For Negatively Charged Decay Vertex Daughter Tracks In TPC", kTH1F, {axisNsigmaE});
    histosConversionPhoton.add("DCAv0Daughters", "DCA Between Decay Vertex Daughter Tracks", kTH1F, {{100, 0.0f, v0setting_dcav0dau, "DCA [cm]"}});
    histosConversionPhoton.add("DCAv0Pos", "DCA Between Decay Vertex And Positively Charged Daughter Track", kTH1F, {axisDCA});
    histosConversionPhoton.add("DCAv0Neg", "DCA Between Decay Vertex And Negatively Charged Daughter Track", kTH1F, {axisDCA});
    histosConversionPhoton.add("v0cosPA", "Cosine Of Decay Vertex Pointing Angle", kTH1F, {{100, v0setting_cospa, 1.0f, "cos#theta"}});
    histosConversionPhoton.add("DCAv0Radius", "Decay Vertex Transverse Radius", kTH1F, {axisV0Radius});
    histosConversionPhoton.add("photonArmenterosPodolanski", "Armenteros-Podolanski Plot For Photon Candidates", kTH2F, {axisAlpha, axisQt}); 

    histosConversionPhoton.add("RecGenPhotondEvsGenE", "Energy Difference Between Reconstructed And Generated Photon vs. Energy", kTH2F, {{100, -0.5f, 0.5f, "E_{rec} - E_{gen} [GeV]"}, {100, 0.0f, maxPhotonPt, "E_{gen} [GeV]"}}); 
    histosConversionPhoton.add("RecGenPhotondPtvsGenPt", "Transverse Momentum Difference Between Reconstructed And Generated Photon vs. Transverse Momentum", kTH2F, {{100, -0.5f, 0.5f, "p_{T, rec} - p_{T, gen} [GeV]"}, {100, 0.0f, maxPhotonPt, "p_{T, gen} [GeV]"}}); 
    histosConversionPhoton.add("MismatchPhoton", "Transverse Momentum Of Mismatched Photons NOT From #Sigma^{0} Hyperons", kTH1F, {{100, 0.0f, maxPhotonPt, "p_{T} [GeV]"}});
    histosConversionPhoton.add("FalsePhoton", "Transverse Momentum Of Mismatched Decay Verteces NOT From Photons", kTH1F, {{100, 0.0f, maxPhotonPt, "p_{T} [GeV]"}});
  }

  void process(filteredCollision const& collision, filteredV0s const& v0s, filteredPhotonDaughterTracks const&, McParticles const&) {
    for (auto const& v0: v0s) {
      histosConversionPhoton.fill(HIST("photonCutEffects"), 0);
      if (v0.v0cosPA() < v0setting_cospa) continue;
      histosConversionPhoton.fill(HIST("photonCutEffects"), 1);
      if (v0.v0radius() < v0setting_radius) continue;
      histosConversionPhoton.fill(HIST("photonCutEffects"), 2);
      if (abs(v0.alpha()) > maxPhotonAlpha) continue;
      histosConversionPhoton.fill(HIST("photonCutEffects"), 3);
      if (v0.qtarm() > maxPhotonQt) continue;
      histosConversionPhoton.fill(HIST("photonCutEffects"), 4);

      auto const& posPhotonDaughterTrack = v0.posTrack_as<filteredPhotonDaughterTracks>();
      auto const& negPhotonDaughterTrack = v0.negTrack_as<filteredPhotonDaughterTracks>();

      float const tpcNPosSigmaEl = posPhotonDaughterTrack.tpcNSigmaEl();
      float const tpcNNegSigmaEl = negPhotonDaughterTrack.tpcNSigmaEl();
      if (abs(tpcNPosSigmaEl) > maxTpcNSigmaEl || abs(tpcNNegSigmaEl) > maxTpcNSigmaEl) continue;
      histosConversionPhoton.fill(HIST("photonCutEffects"), 5);

      float const photonPx = posPhotonDaughterTrack.px() + negPhotonDaughterTrack.px();
      float const photonPy = posPhotonDaughterTrack.py() + negPhotonDaughterTrack.py();
      float const photonPz = posPhotonDaughterTrack.pz() + negPhotonDaughterTrack.pz();
      float const electronEnergy = std::sqrt(negPhotonDaughterTrack.p()*negPhotonDaughterTrack.p() + electronMass*electronMass);
      float const positronEnergy = std::sqrt(posPhotonDaughterTrack.p()*posPhotonDaughterTrack.p() + electronMass*electronMass);
      float const photonEnergy = electronEnergy + positronEnergy;
      TLorentzVector const photon(photonPx, photonPy, photonPz, photonEnergy);
      if (photon.Pt() > maxPhotonPt) continue;
      histosConversionPhoton.fill(HIST("photonCutEffects"), 6);
      if (photon.M() > maxPhotonMass || photon.M() < 0.0f) continue;
      histosConversionPhoton.fill(HIST("photonCutEffects"), 7);

      histosConversionPhoton.fill(HIST("posTPCEl"), tpcNPosSigmaEl);
      histosConversionPhoton.fill(HIST("negTPCEl"), tpcNNegSigmaEl);

      histosConversionPhoton.fill(HIST("DCAv0Daughters"), v0.dcaV0daughters());
      histosConversionPhoton.fill(HIST("DCAv0Pos"), v0.dcapostopv());
      histosConversionPhoton.fill(HIST("DCAv0Neg"), v0.dcanegtopv());
      histosConversionPhoton.fill(HIST("v0cosPA"), v0.v0cosPA());
      histosConversionPhoton.fill(HIST("DCAv0Radius"), v0.v0radius());
      histosConversionPhoton.fill(HIST("photonArmenterosPodolanski"), v0.alpha(), v0.qtarm());

      histosConversionPhoton.fill(HIST("photonMass"), photon.M());
      histosConversionPhoton.fill(HIST("photonPt"), photon.Pt());
      histosConversionPhoton.fill(HIST("photonEta"), photon.Eta());

      AddConversionPhoton(photon.Px(), photon.Py(), photon.Pz(), photon.E());

      if (v0.has_mcParticle()) {
        auto const& v0mcParticle = v0.mcParticle();
        // check that the V0 comes from a photon
        if (v0mcParticle.pdgCode() == 22) {
          if (v0mcParticle.has_mothers()) {
            for (auto const& photonMother: v0mcParticle.mothers_as<McParticles>()) {
              // check that photon comes from Sigma^0 hyperon:
              if (abs(photonMother.pdgCode()) == 3212) {
                histosConversionPhoton.fill(HIST("RecGenPhotondEvsGenE"), photonEnergy - v0mcParticle.e(), v0mcParticle.e());
                histosConversionPhoton.fill(HIST("RecGenPhotondPtvsGenPt"), photon.Pt() - v0mcParticle.pt(), v0mcParticle.pt());
              } else {
                histosConversionPhoton.fill(HIST("MismatchPhoton"), v0mcParticle.pt());    
              }
            }
          }
        } else {
          histosConversionPhoton.fill(HIST("FalsePhoton"), v0mcParticle.pt());    
        }
      }
    }   
  }
};

struct ProcessLambdaHyperons {
  Produces<LambdaHyperon> AddLambdaHyperon;

  float chargedPionMass = o2::constants::physics::MassPionCharged;
  float protonMass = o2::constants::physics::MassProton;

  Configurable<float> zVertexCut{"zVertexCut", 10.0f, "Maximum Primary Vertex Z coordinate [cm]"};
  Filter zVertexFilter = (nabs(collision::posZ) < zVertexCut);
  Filter zVertexErrorFilter = (collision::posZ != 0.0f);
  Filter eventSelectionFilter = (evsel::sel8 == true);
  using filteredCollision = Filtered<Join<Collisions, EvSels>>::iterator;
    
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsMass{"nBinsMass", 100, "N bins in invariant mass histo"};
  Configurable<float> etaCut{"etaCut", 1.2f, "Maximum Pseudorapidity"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1.0f, "DCA V0 Daughters [cm]"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06f, "DCA Pos To PV [cm]"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06f, "DCA Neg To PV [cm]"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.995, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.0f, "V0 Radius [cm]"};

  Configurable<float> maxLambdaAlpha{"maxLambdaAlpha", 0.9f, "Maximum Lambda Decay Asymmetry"};
  Configurable<float> minLambdaAlpha{"minLambdaAlpha", 0.4f, "Minimum Lambda Decay Asymmetry"};
  Configurable<float> maxLambdaQt{"maxLambdaQt", 0.13f, "Maximum Lambda Decay q_{T} [GeV]"};
  Configurable<float> minLambdaQt{"minLambdaQt", 0.02f, "Minimum Lambda Decay  [GeV]"};

  // filter can only be applied to static columns
  Filter preV0Filter = (nabs(v0data::dcapostopv) > v0setting_dcapostopv && nabs(v0data::dcanegtopv) > v0setting_dcanegtopv && v0data::dcaV0daughters < v0setting_dcav0dau);
  using filteredV0s = Filtered<Join<V0Datas, McV0Labels>>;

  // Tracks-related:
  Configurable<float> minLambdaMass{"minLambdaMass", 1.10f, "Minimum Lambda Hyperon Invariant Mass [GeV]"};
  Configurable<float> maxLambdaMass{"maxLambdaMass", 1.13f, "Maximum Lambda Hyperon Invariant Mass [GeV]"};
  Configurable<float> maxLambdaPt{"maxLambdaPt", 5.0f, "Maximum Lambda Hyperon Transverse Momentum [GeV]"};

  Configurable<float> maxTpcNSigmaPr{"maxTpcNSigmaPr", 3.0f, "Maximum N_{#sigma_{p}} from TPC signal"};
  Configurable<float> maxTpcNSigmaPi{"maxTpcNSigmaPi", 3.0f, "Maximum N_{#sigma_{#pi}} from TPC signal"};

  Configurable<float> maxTrackDCA{"maxTrackDCA", 1.0f, "Maximum Track DCA [cm]"};

  Filter etaFilter = (nabs(track::eta) < etaCut);
  Filter dcaFilter = (nabs(track::dcaXY) < maxTrackDCA);

  using LambdaDaughterTracks = Join<Tracks, TracksDCA, pidTPCPi, pidTPCPr, McTrackLabels>;
  using filteredLambdaDaughterTracks = Filtered<LambdaDaughterTracks>;

  HistogramRegistry histosLambda{"histosLambda", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&) {
    const AxisSpec axisLambdaMass{nBinsMass, minLambdaMass, maxLambdaMass, "M_{p#pi} [GeV]"};
    const AxisSpec axisLambdaPt{nBinsPt, 0.0f, maxLambdaPt, "p_{T} [GeV]"};
    const AxisSpec axisEta{150, -etaCut, etaCut, "#eta"};
    const AxisSpec axisDCA{100, -5.0f, 5.0f, "DCA [cm]"};
    const AxisSpec axisV0Radius{100, 0.0f, 50.0f, "V0 Radius [cm]"};
    const AxisSpec axisAlpha{100, -1.0f, 1.0f, "#alpha"};
    const AxisSpec axisQt{100, 0.0f, 0.25f, "q_{T} [GeV]"};
    const AxisSpec axisNsigmaPr{100, -maxTpcNSigmaPr, maxTpcNSigmaPr, "N_{#sigma_{p}}"};
    const AxisSpec axisNsigmaPi{100, -maxTpcNSigmaPi, maxTpcNSigmaPi, "N_{#sigma_{#pi}}"};

    histosLambda.add("LambdaCutEffects", "Cuts Effect On Number Of Kept #Lambda Hyperon Candidates", kTH1D, {{10, 0, 10, "Cut"}});

    histosLambda.add("LambdaEta", "Pseudorapidity Of #Lambda Hyperon Candidates", kTH1F, {axisEta});
    histosLambda.add("LambdaMass", "Invariant Mass Of #Lambda Hyperon Candidates", kTH1F, {axisLambdaMass});
    histosLambda.add("LambdaPt", "Transverse Momentum Of #Lambda Hyperon Candidates", kTH1F, {axisLambdaPt});

    histosLambda.add("posTPCPr", "N_{#sigma, p} For Positively Charged Decay Vertex Daughter Tracks In TPC", kTH1F, {axisNsigmaPr});
    histosLambda.add("negTPCPr", "N_{#sigma, p} For Negatively Charged Decay Vertex Daughter Tracks In TPC", kTH1F, {axisNsigmaPr});     
    histosLambda.add("posTPCPi", "N_{#sigma, #pi} For Positively Charged Decay Vertex Daughter Tracks In TPC", kTH1F, {axisNsigmaPi});
    histosLambda.add("negTPCPi", "N_{#sigma, #pi} For Negatively Charged Decay Vertex Daughter Tracks In TPC", kTH1F, {axisNsigmaPi}); 
    histosLambda.add("DCAv0Daughters", "DCA Between Decay Vertex Daughter Tracks", kTH1F, {{100, 0.0f, v0setting_dcav0dau, "DCA [cm]"}});
    histosLambda.add("DCAv0Pos", "DCA Between Decay Vertex And Positively Charged Daughter Track", kTH1F, {axisDCA});
    histosLambda.add("DCAv0Neg", "DCA Between Decay Vertex And Negatively Charged Daughter Track", kTH1F, {axisDCA});
    histosLambda.add("v0cosPA", "Cosine Of Decay Vertex Pointing Angle", kTH1F, {{100, v0setting_cospa, 1.0f, "cos#theta"}});
    histosLambda.add("DCAv0Radius", "Decay Vertex Transverse Radius", kTH1F, {axisV0Radius});
    histosLambda.add("LambdaArmenterosPodolanski", "Armenteros-Podolanski Plot For #Lambda Hyperon Candidates", kTH2F, {axisAlpha, axisQt}); 

    histosLambda.add("RecGenLambdadEvsGenE", "Energy Difference Between Reconstructed And Generated #Lambda Hyperon vs. Energy", kTH2F, {{100, -0.5f, 0.5f, "E_{rec} - E_{gen} [GeV]"}, {100, 0.0f, maxLambdaPt, "E_{gen} [GeV]"}}); 
    histosLambda.add("RecGenLambdadPtvsGenPt", "Transverse Momentum Difference Between Reconstructed And Generated #Lambda Hyperon vs. Transverse Momentum", kTH2F, {{100, -0.5f, 0.5f, "p_{T, rec} - p_{T, gen} [GeV]"}, {100, 0.0f, maxLambdaPt, "p_{T, gen} [GeV]"}}); 
    histosLambda.add("MismatchLambda", "Transverse Momentum Of Mismatched #Lambda Hyperons NOT From #Sigma^{0} Hyperons", kTH1F, {{100, 0.0f, maxLambdaPt, "p_{T} [GeV]"}});
    histosLambda.add("FalseLambda", "Transverse Momentum Of Mismatched Decay Verteces NOT From #Lambda Hyperons", kTH1F, {{100, 0.0f, maxLambdaPt, "p_{T} [GeV]"}});
  }

  void process(filteredCollision const& collision, filteredV0s const& v0s, filteredLambdaDaughterTracks const&, McParticles const&) {
    for (auto const& v0: v0s) {
      histosLambda.fill(HIST("LambdaCutEffects"), 0);
      if (v0.v0cosPA() < v0setting_cospa) continue;
      histosLambda.fill(HIST("LambdaCutEffects"), 1);
      if (v0.v0radius() < v0setting_radius) continue;
      histosLambda.fill(HIST("LambdaCutEffects"), 2);
      if (abs(v0.alpha()) > maxLambdaAlpha || abs(v0.alpha()) < minLambdaAlpha) continue;
      histosLambda.fill(HIST("LambdaCutEffects"), 3);
      if (v0.qtarm() > maxLambdaQt || v0.qtarm() < minLambdaQt) continue;
      histosLambda.fill(HIST("LambdaCutEffects"), 4);

      auto const& posLambdaDaughterTrack = v0.posTrack_as<filteredLambdaDaughterTracks>();
      auto const& negLambdaDaughterTrack = v0.negTrack_as<filteredLambdaDaughterTracks>();

      float const tpcNPosSigmaPr = posLambdaDaughterTrack.tpcNSigmaPr();
      float const tpcNPosSigmaPi = posLambdaDaughterTrack.tpcNSigmaPi();

      float const tpcNNegSigmaPr = negLambdaDaughterTrack.tpcNSigmaPr();
      float const tpcNNegSigmaPi = negLambdaDaughterTrack.tpcNSigmaPi();

      bool posTraskIsProton = abs(tpcNPosSigmaPr) < maxTpcNSigmaPr;
      bool posTraskIsPion = abs(tpcNPosSigmaPi) < maxTpcNSigmaPi;

      bool negTraskIsProton = abs(tpcNNegSigmaPr) < maxTpcNSigmaPr;
      bool negTraskIsPion = abs(tpcNNegSigmaPi) < maxTpcNSigmaPi;

      bool isLambda = false;
      bool isAntiLambda = false;

      if (posTraskIsProton && negTraskIsPion) isLambda = true;
      if (negTraskIsProton && posTraskIsPion) isAntiLambda = true;

      if (!isLambda && !isAntiLambda) continue;
      histosLambda.fill(HIST("LambdaCutEffects"), 5);

      float const LambdaPx = posLambdaDaughterTrack.px() + negLambdaDaughterTrack.px();
      float const LambdaPy = posLambdaDaughterTrack.py() + negLambdaDaughterTrack.py();
      float const LambdaPz = posLambdaDaughterTrack.pz() + negLambdaDaughterTrack.pz();
      float protonEnergy = 0;
      float pionEnergy = 0;
      if (isLambda) {
        protonEnergy = std::sqrt(posLambdaDaughterTrack.p()*posLambdaDaughterTrack.p() + protonMass*protonMass);
        pionEnergy = std::sqrt(negLambdaDaughterTrack.p()*negLambdaDaughterTrack.p() + chargedPionMass*chargedPionMass);
      } else if (isAntiLambda) {
        protonEnergy = std::sqrt(negLambdaDaughterTrack.p()*negLambdaDaughterTrack.p() + protonMass*protonMass);
        pionEnergy = std::sqrt(posLambdaDaughterTrack.p()*posLambdaDaughterTrack.p() + chargedPionMass*chargedPionMass);
      }
      float const LambdaEnergy = protonEnergy + pionEnergy;
      TLorentzVector const Lambda(LambdaPx, LambdaPy, LambdaPz, LambdaEnergy);
      if (Lambda.Pt() > maxLambdaPt) continue;
      histosLambda.fill(HIST("LambdaCutEffects"), 6);
      if (Lambda.M() > maxLambdaMass || Lambda.M() < minLambdaMass) continue;
      histosLambda.fill(HIST("LambdaCutEffects"), 7);

      if (isLambda) {
        histosLambda.fill(HIST("posTPCPr"), tpcNPosSigmaPr);
        histosLambda.fill(HIST("negTPCPi"), tpcNNegSigmaPi);
      } else if (isAntiLambda) {
        histosLambda.fill(HIST("posTPCPi"), tpcNPosSigmaPi);
        histosLambda.fill(HIST("negTPCPr"), tpcNNegSigmaPr);
      }

      histosLambda.fill(HIST("DCAv0Daughters"), v0.dcaV0daughters());
      histosLambda.fill(HIST("DCAv0Pos"), v0.dcapostopv());
      histosLambda.fill(HIST("DCAv0Neg"), v0.dcanegtopv());
      histosLambda.fill(HIST("v0cosPA"), v0.v0cosPA());
      histosLambda.fill(HIST("DCAv0Radius"), v0.v0radius());
      histosLambda.fill(HIST("LambdaArmenterosPodolanski"), v0.alpha(), v0.qtarm());

      histosLambda.fill(HIST("LambdaMass"), Lambda.M());
      histosLambda.fill(HIST("LambdaPt"), Lambda.Pt());
      histosLambda.fill(HIST("LambdaEta"), Lambda.Eta());

      AddLambdaHyperon(Lambda.Px(), Lambda.Py(), Lambda.Pz(), Lambda.E());

      if (v0.has_mcParticle()) {
        auto const& v0mcParticle = v0.mcParticle();
        // check that the V0 comes from a lambda hyperon
        if (abs(v0mcParticle.pdgCode()) == 3122) {
          if (v0mcParticle.has_mothers()) {
            for (auto const& LambdaMother: v0mcParticle.mothers_as<McParticles>()) {
              // check that Lambda hyperon comes from Sigma^0 hyperon:
              if (abs(LambdaMother.pdgCode()) == 3212) { 
                histosLambda.fill(HIST("RecGenLambdadEvsGenE"), LambdaEnergy - v0mcParticle.e(), v0mcParticle.e());
                histosLambda.fill(HIST("RecGenLambdadPtvsGenPt"), Lambda.Pt() - v0mcParticle.pt(), v0mcParticle.pt());
              } else {
                histosLambda.fill(HIST("MismatchLambda"), v0mcParticle.pt());       
              }
            }              
          }
        } else {
          histosLambda.fill(HIST("FalseLambda"), v0mcParticle.pt());       
        }
      }
    }   
  }
};

struct ReconstructSigma0viaPCM {
  Configurable<float> zVertexCut{"zVertexCut", 10.0f, "Maximum Primary Vertex Z coordinate [cm]"};
  Filter zVertexFilter = (nabs(collision::posZ) < zVertexCut);
  Filter zVertexErrorFilter = (collision::posZ != 0.0f);
  Filter eventSelectionFilter = (evsel::sel8 == true);
  using filteredCollision = Filtered<Join<Collisions, EvSels>>::iterator;

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsMass{"nBinsMass", 100, "N bins in invariant mass histo"};
  Configurable<float> minSigma0Mass{"minSigma0Mass", 1.1f, "Maximum Sigma^0 Invariant Mass [GeV]"};
  Configurable<float> maxSigma0Mass{"maxSigma0Mass", 1.3f, "Maximum Sigma^0 Invariant Mass [GeV]"};
  Configurable<float> maxSigma0Pt{"maxSigma0Pt", 5.0f, "Maximum Sigma^0 Hyperon Transverse Momentum [GeV]"};

  HistogramRegistry histosSigma0{"histosSigma0", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&) {
    const AxisSpec axisSigma0Mass{nBinsMass, minSigma0Mass, maxSigma0Mass, "M_{#Lambda#gamma_{PCM}} [GeV]"};
    const AxisSpec axisSigma0Pt{nBinsPt, 0.0f, maxSigma0Pt, "p_{T, candidate} [GeV]"};

    histosSigma0.add("Sigma0PCM", "#Sigma^{0} Candidates From Conversion Photons", kTH2F, {axisSigma0Mass, axisSigma0Pt});
  }

  void process(filteredCollision const&, ConversionPhoton const& photons, LambdaHyperon const& Lambdas) {
    for (auto const& [Lambda, photon]: combinations(CombinationsFullIndexPolicy(Lambdas, photons))) {
      TLorentzVector const LambdaLorentzVector(Lambda.px(), Lambda.py(), Lambda.pz(), Lambda.e());
      TLorentzVector const photonLorentzVector(photon.px(), photon.py(), photon.pz(), photon.e());
      TLorentzVector const Sigma0LorentzVector = LambdaLorentzVector + photonLorentzVector;
      float const Sigma0Mass = Sigma0LorentzVector.M();
      float const Sigma0Pt = Sigma0LorentzVector.Pt();
      if (Sigma0Mass < maxSigma0Mass && Sigma0Mass > minSigma0Mass) {
        histosSigma0.fill(HIST("Sigma0PCM"), Sigma0Mass, Sigma0Pt);           
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{
    adaptAnalysisTask<ProcessCollisions>(cfgc),
    adaptAnalysisTask<ProcessGeneratedEvents>(cfgc),
    adaptAnalysisTask<ProcessConversionPhotons>(cfgc),
    adaptAnalysisTask<ProcessLambdaHyperons>(cfgc),
    adaptAnalysisTask<ReconstructSigma0viaPCM>(cfgc)
  };
}