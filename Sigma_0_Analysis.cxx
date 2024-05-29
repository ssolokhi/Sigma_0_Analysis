// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "LICENSE".
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
  namespace mcconversionphoton {
    DECLARE_SOA_COLUMN(GenPx, genpx, float);
    DECLARE_SOA_COLUMN(RecPx, recpx, float);
    DECLARE_SOA_COLUMN(GenPy, genpy, float);
    DECLARE_SOA_COLUMN(RecPy, recpy, float);
    DECLARE_SOA_COLUMN(GenPz, genpz, float);
    DECLARE_SOA_COLUMN(RecPz, recpz, float);
    DECLARE_SOA_COLUMN(GenE, gene, float);
    DECLARE_SOA_COLUMN(RecE, rece, float);
  }

  DECLARE_SOA_TABLE(mcConversionPhoton, "AOD", "MCPCMPHOTON", mcconversionphoton::GenPx, mcconversionphoton::RecPx,
  mcconversionphoton::GenPy, mcconversionphoton::RecPy, mcconversionphoton::GenPz, mcconversionphoton::RecPz,
  mcconversionphoton::GenE, mcconversionphoton::RecE);
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

namespace o2::aod {
  namespace mclambdahyperon {
    DECLARE_SOA_COLUMN(GenPx, genpx, float);
    DECLARE_SOA_COLUMN(RecPx, recpx, float);
    DECLARE_SOA_COLUMN(GenPy, genpy, float);
    DECLARE_SOA_COLUMN(RecPy, recpy, float);
    DECLARE_SOA_COLUMN(GenPz, genpz, float);
    DECLARE_SOA_COLUMN(RecPz, recpz, float);
    DECLARE_SOA_COLUMN(GenE, gene, float);
    DECLARE_SOA_COLUMN(RecE, rece, float);
  }

  DECLARE_SOA_TABLE(mcLambdaHyperon, "AOD", "MCLAMBDAHYPERON", mclambdahyperon::GenPx, mclambdahyperon::RecPx,
  mclambdahyperon::GenPy, mclambdahyperon::RecPy, mclambdahyperon::GenPz, mclambdahyperon::RecPz,
  mclambdahyperon::GenE, mclambdahyperon::RecE);
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
    const AxisSpec axisDCADaughters{100, 0.0f, 1.0f, "DCA [cm]"};
    const AxisSpec axisDCA{100, -5.0f, 5.0f, "DCA [cm]"};
    const AxisSpec axisV0Radius{100, 0.0f, 50.0f, "V0 Radius [cm]"};
    const AxisSpec axisAlpha{100, -1.0f, 1.0f, "#alpha"};
    const AxisSpec axisQt{100, 0.0f, 0.25f, "q_{T} [GeV]"};
    const AxisSpec axisCosTheta{100, 0.99f, 1.0f, "cos#theta"};

    histosMC.add("mcGenEta", "Pseudorapidity Of Generated Tracks (MC)", kTH1F, {axisEta});
    histosMC.add("mcGenPt", "Transverse Momentum Of Generated Tracks (MC)", kTH1F, {{100, 0.0f, 10.0f, "p_{T} [GeV]"}});
    histosMC.add("mcGenPhotonPt", "Transverse Momentum Of Generated Photons (MC)", kTH1F, {axisPhotonPt});    
    histosMC.add("mcGenLambdaPt", "Transverse Momentum Of Generated #Lambda Hyperons (MC)", kTH1F, {axisLambdaPt});
    histosMC.add("mcGenSigma0Pt", "Transverse Momentum Of Generated #Sigma^{0} Hyperons (MC)", kTH1F, {axisSigma0Pt});
    histosMC.add("mcGenPhotonFromSigma0Pt", "Transverse Momentum Of Generated Photons From #Sigma^{0} Hyperons (MC)", kTH1F, {axisPhotonPt});
    histosMC.add("mcGenLambdaFromSigma0Pt", "Transverse Momentum Of Generated #Lambda Hyperons From #Sigma^{0} Hyperons (MC)", kTH1F, {axisLambdaPt});
  
    histosMC.add("mcPhotonDCAv0Daughters", "DCA Between Photon Vertex Daughter Tracks (MC)", kTH1F, {axisDCADaughters});
    histosMC.add("mcLambdaDCAv0Daughters", "DCA Between #Lambda Hyperon Vertex Daughter Tracks (MC)", kTH1F, {axisDCADaughters});
    histosMC.add("mcOtherDCAv0Daughters", "DCA Between Background Vertex Daughter Tracks (MC)", kTH1F, {axisDCADaughters});

    histosMC.add("mcPhotonDCAv0Pos", "DCA Between Photon Vertex And Positively Charged Daughter Track (MC)", kTH1F, {axisDCA});
    histosMC.add("mcPhotonDCAv0Neg", "DCA Between Photon Vertex And Negatively Charged Daughter Track (MC)", kTH1F, {axisDCA});
    histosMC.add("mcLambdaDCAv0Pos", "DCA Between #Lambda Hyperon And Positively Charged Daughter Track (MC)", kTH1F, {axisDCA});
    histosMC.add("mcLambdaDCAv0Neg", "DCA Between #Lambda Hyperon And Negatively Charged Daughter Track (MC)", kTH1F, {axisDCA});
    histosMC.add("mcOtherDCAv0Pos", "DCA Between Background Vertex And Positively Charged Daughter Track (MC)", kTH1F, {axisDCA});
    histosMC.add("mcOtherDCAv0Neg", "DCA Between Background Vertex And Negatively Charged Daughter Track (MC)", kTH1F, {axisDCA});

    histosMC.add("mcPhotonV0cosPA", "Cosine Of Photon Vertex Pointing Angle (MC)", kTH1F, {axisCosTheta});
    histosMC.add("mcLambdaV0cosPA", "Cosine Of #Lambda Hyperon Vertex Pointing Angle (MC)", kTH1F, {{axisCosTheta}});
    histosMC.add("mcOtherV0cosPA", "Cosine Of Background Vertex Pointing Angle (MC)", kTH1F, {axisCosTheta});

    histosMC.add("mcPhotonV0Radius", "Photon Vertex Transverse Radius (MC)", kTH1F, {axisV0Radius});
    histosMC.add("mcLambdaV0Radius", "#Lambda Hyperon Vertex Transverse Radius (MC)", kTH1F, {axisV0Radius});
    histosMC.add("mcOtherV0Radius", "Background Vertex Transverse Radius (MC)", kTH1F, {axisV0Radius});

    histosMC.add("mcPhotonArmenterosPodolanski", "Armenteros-Podolanski Plot For Photons (MC)", kTH2F, {axisAlpha, axisQt}); 
    histosMC.add("mcLambdaArmenterosPodolanski", "Armenteros-Podolanski Plot For #Lambda Hyperons (MC)", kTH2F, {axisAlpha, axisQt}); 
    histosMC.add("mcOtherArmenterosPodolanski", "Armenteros-Podolanski Plot For Background Candidates (MC)", kTH2F, {axisAlpha, axisQt}); 
  }

  void process(McCollision const&, McParticles const& mcparticles, Join<V0Datas, McV0Labels> const& v0s) {
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

    for (auto const& v0: v0s) {
      if (v0.has_mcParticle()) {
        auto const& v0mcParticle = v0.mcParticle();
        if (abs(v0mcParticle.pdgCode()) == 3122) {
          histosMC.fill(HIST("mcLambdaDCAv0Daughters"), v0.dcaV0daughters());
          histosMC.fill(HIST("mcLambdaDCAv0Pos"), v0.dcapostopv());
          histosMC.fill(HIST("mcLambdaDCAv0Neg"), v0.dcanegtopv());
          histosMC.fill(HIST("mcLambdaV0cosPA"), v0.v0cosPA());
          histosMC.fill(HIST("mcLambdaV0Radius"), v0.v0radius());
          histosMC.fill(HIST("mcLambdaArmenterosPodolanski"), v0.alpha(), v0.qtarm());
        } else if (v0mcParticle.pdgCode() == 22) {
          histosMC.fill(HIST("mcPhotonDCAv0Daughters"), v0.dcaV0daughters());
          histosMC.fill(HIST("mcPhotonDCAv0Pos"), v0.dcapostopv());
          histosMC.fill(HIST("mcPhotonDCAv0Neg"), v0.dcanegtopv());
          histosMC.fill(HIST("mcPhotonV0cosPA"), v0.v0cosPA());
          histosMC.fill(HIST("mcPhotonV0Radius"), v0.v0radius());
          histosMC.fill(HIST("mcPhotonArmenterosPodolanski"), v0.alpha(), v0.qtarm());
        } else {
          histosMC.fill(HIST("mcOtherDCAv0Daughters"), v0.dcaV0daughters());
          histosMC.fill(HIST("mcOtherDCAv0Pos"), v0.dcapostopv());
          histosMC.fill(HIST("mcOtherDCAv0Neg"), v0.dcanegtopv());
          histosMC.fill(HIST("mcOtherV0cosPA"), v0.v0cosPA());
          histosMC.fill(HIST("mcOtherV0Radius"), v0.v0radius());
          histosMC.fill(HIST("mcOtherArmenterosPodolanski"), v0.alpha(), v0.qtarm());
        }
      }
    }
  }
};

struct ProcessConversionPhotons {
  Produces<ConversionPhoton> AddConversionPhoton;
  Produces<mcConversionPhoton> AddMCconversionPhoton;

  float electronMass = o2::constants::physics::MassElectron;

  Configurable<float> zVertexCut{"zVertexCut", 10.0f, "Maximum Primary Vertex Z coordinate [cm]"};
  Filter zVertexFilter = (nabs(collision::posZ) < zVertexCut);
  Filter zVertexErrorFilter = (collision::posZ != 0.0f);
  //Filter eventSelectionFilter = (evsel::sel8 == true);
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

  Configurable<float> maxPhotonMass{"maxPhotonMass", 1.0f, "Maximum Electron-Positron Invariant Mass [GeV]"};
  Configurable<float> maxPhotonPt{"maxPhotonPt", 10.0f, "Maximum Photon Transverse Momentum [GeV]"};
  Configurable<float> maxTpcNSigmaEl{"maxTpcNSigmaEl", 3.0f, "Maximum N_{#sigma_{e}} from TPC signal"};
  Configurable<float> maxTrackDCA{"maxTrackDCA", 1.0f, "Maximum Track DCA [cm]"};
  Configurable<float> maxPhotonAlpha{"maxPhotonAlpha", 0.8f, "Maximum Photon Decay Asymmetry"};
  Configurable<float> maxPhotonQt{"maxPhotonQt", 0.1f, "Maximum Photon Decay q_{T} [GeV]"};

  Filter etaFilter = (nabs(track::eta) < etaCut);
  Filter dcaFilter = (nabs(track::dcaXY) < maxTrackDCA);

  using photonDaughterTracks = Join<Tracks, TracksDCA, pidTPCEl, McTrackLabels>;
  using filteredPhotonDaughterTracks = Filtered<photonDaughterTracks>;

  HistogramRegistry histosConversionPhoton{"histosConversionPhoton", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&) {
    printf("Init photon");
    const AxisSpec axisEta{100, -etaCut, etaCut, "#eta"};
    const AxisSpec axisPhotonMass{nBinsMass, 0.0f, maxPhotonMass, "M_{e^{+}e^{-}} [GeV]"};
    const AxisSpec axisPhotonPt{nBinsPt, 0.0f, maxPhotonPt, "p_{T} [GeV]"};
    const AxisSpec axisNsigmaE{100, -maxTpcNSigmaEl, maxTpcNSigmaEl, "N_{#sigma_{e}}"};
    const AxisSpec axisDCA{100, -5.0f, 5.0f, "DCA [cm]"};
    const AxisSpec axisV0Radius{100, 0.0f, 100.0f, "V0 Radius [cm]"};
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
    histosConversionPhoton.add("MismatchedPhoton", "Transverse Momentum Of Mismatched Photons NOT From #Sigma^{0} Hyperons", kTH1F, {{100, 0.0f, maxPhotonPt, "p_{T} [GeV]"}});
    histosConversionPhoton.add("FalsePhoton", "Transverse Momentum Of Mismatched Decay Verteces NOT From Photons", kTH1F, {{100, 0.0f, maxPhotonPt, "p_{T} [GeV]"}});
  }

  void process(filteredCollision const& collision, filteredV0s const& v0s, filteredPhotonDaughterTracks const&, McParticles const&) {
    for (auto const& v0: v0s) {
      printf("Photon\n");
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
            auto const& photonMother = v0mcParticle.mothers_first_as<McParticles>();
            // check that photon comes from Sigma^0 hyperon:
            if (abs(photonMother.pdgCode()) == 3212) {
              AddMCconversionPhoton(v0mcParticle.px(), photon.Px(), v0mcParticle.py(), photon.Py(),
              v0mcParticle.pz(), photon.Py(), v0mcParticle.e(), photon.E());
              histosConversionPhoton.fill(HIST("RecGenPhotondEvsGenE"), photonEnergy - v0mcParticle.e(), v0mcParticle.e());
              histosConversionPhoton.fill(HIST("RecGenPhotondPtvsGenPt"), photon.Pt() - v0mcParticle.pt(), v0mcParticle.pt());
            } else {
              histosConversionPhoton.fill(HIST("MismatchedPhoton"), v0mcParticle.pt());    
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
  Produces<mcLambdaHyperon> AddMCLambdaHyperon;

  float chargedPionMass = o2::constants::physics::MassPionCharged;
  float protonMass = o2::constants::physics::MassProton;

  Configurable<float> zVertexCut{"zVertexCut", 10.0f, "Maximum Primary Vertex Z coordinate [cm]"};
  Filter zVertexFilter = (nabs(collision::posZ) < zVertexCut);
  Filter zVertexErrorFilter = (collision::posZ != 0.0f);
  //Filter eventSelectionFilter = (evsel::sel8 == true);
  using filteredCollision = Filtered<Join<Collisions, EvSels>>::iterator;
    
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsMass{"nBinsMass", 100, "N bins in invariant mass histo"};
  Configurable<float> etaCut{"etaCut", 1.2f, "Maximum Pseudorapidity"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 2.0f, "DCA V0 Daughters [cm]"};
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
  Configurable<float> maxLambdaPt{"maxLambdaPt", 10.0f, "Maximum Lambda Hyperon Transverse Momentum [GeV]"};

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
    histosLambda.add("MismatchedLambda", "Transverse Momentum Of Mismatched #Lambda Hyperons NOT From #Sigma^{0} Hyperons", kTH1F, {{100, 0.0f, maxLambdaPt, "p_{T} [GeV]"}});
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
            auto const& LambdaMother = v0mcParticle.mothers_first_as<McParticles>();
            // check that Lambda hyperon comes from Sigma^0 hyperon:
            if (abs(LambdaMother.pdgCode()) == 3212) { 
              AddMCLambdaHyperon(v0mcParticle.px(), Lambda.Px(), v0mcParticle.py(), Lambda.Py(),
              v0mcParticle.pz(), Lambda.Py(), v0mcParticle.e(), Lambda.E());
              histosLambda.fill(HIST("RecGenLambdadEvsGenE"), LambdaEnergy - v0mcParticle.e(), v0mcParticle.e());
              histosLambda.fill(HIST("RecGenLambdadPtvsGenPt"), Lambda.Pt() - v0mcParticle.pt(), v0mcParticle.pt());
            } else {
              histosLambda.fill(HIST("MismatchedLambda"), v0mcParticle.pt());       
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
  //Filter eventSelectionFilter = (evsel::sel8 == true);
  using filteredCollision = Filtered<Join<Collisions, EvSels>>::iterator;

  Configurable<int> nBinsPt{"nBinsPt", 120, "N bins in pT histo"};
  Configurable<int> nBinsMass{"nBinsMass", 150, "N bins in invariant mass histo"};
  Configurable<float> minSigma0Mass{"minSigma0Mass", 1.12f, "Maximum Sigma^0 Invariant Mass [GeV]"};
  Configurable<float> maxSigma0Mass{"maxSigma0Mass", 1.32f, "Maximum Sigma^0 Invariant Mass [GeV]"};
  Configurable<float> maxSigma0Pt{"maxSigma0Pt", 12.0f, "Maximum Sigma^0 Hyperon Transverse Momentum [GeV]"};

  HistogramRegistry histosSigma0{"histosSigma0", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&) {
    const AxisSpec axisSigma0Mass{nBinsMass, minSigma0Mass, maxSigma0Mass, "M_{#Lambda#gamma_{PCM}} [GeV]"};
    const AxisSpec axisSigma0Pt{nBinsPt, 0.0f, maxSigma0Pt, "p_{T, candidate} [GeV]"};
    const AxisSpec axisAlpha{100, -1.0f, 1.0f, "#alpha"};
    const AxisSpec axisQt{100, 0.0f, 0.25f, "q_{T} [GeV]"};

    histosSigma0.add("Sigma0PCM", "#Sigma^{0} Candidates From Conversion Photons", kTH2F, {axisSigma0Mass, axisSigma0Pt});
    histosSigma0.add("Sigma0PCMArmenterosPodolanski", "Armenteros-Podolanski Plot for #Sigma^{0} Candidates (PCM)", kTH2F, {axisAlpha, axisQt});
  }

  void process(filteredCollision const&, ConversionPhoton const& photons, LambdaHyperon const& Lambdas) {
    for (auto const& [Lambda, photon]: combinations(CombinationsFullIndexPolicy(Lambdas, photons))) {
      TLorentzVector const LambdaLorentzVector(Lambda.px(), Lambda.py(), Lambda.pz(), Lambda.e());
      TLorentzVector const photonLorentzVector(photon.px(), photon.py(), photon.pz(), photon.e());
      TLorentzVector const Sigma0LorentzVector = LambdaLorentzVector + photonLorentzVector;
      float const Sigma0Mass = Sigma0LorentzVector.M();
      float const Sigma0Pt = Sigma0LorentzVector.Pt();

      TVector3 LambdaMomentum(Lambda.px(), Lambda.py(), Lambda.pz());
      TVector3 photonMomentum(photon.px(), photon.py(), photon.pz());
      TVector3 Sigma0Momentum(Sigma0LorentzVector.Px(), Sigma0LorentzVector.Py(), Sigma0LorentzVector.Pz());
      float pLongPos = LambdaMomentum.Dot(Sigma0Momentum)/Sigma0Momentum.Mag();
      float pLongNeg = photonMomentum.Dot(Sigma0Momentum)/Sigma0Momentum.Mag();
      float alpha = (pLongPos - pLongNeg)/(pLongPos + pLongNeg);
      float qT = photonMomentum.Perp(Sigma0Momentum);

      if (Sigma0Mass < maxSigma0Mass && Sigma0Mass > minSigma0Mass) {
        histosSigma0.fill(HIST("Sigma0PCM"), Sigma0Mass, Sigma0Pt);
        histosSigma0.fill(HIST("Sigma0PCMArmenterosPodolanski"), alpha, qT);
      }
    }
  }
};

struct ReconstructNeutralPionsViaPCM {
  Configurable<float> zVertexCut{"zVertexCut", 10.0f, "Maximum Primary Vertex Z coordinate [cm]"};
  Filter zVertexFilter = (nabs(collision::posZ) < zVertexCut);
  Filter zVertexErrorFilter = (collision::posZ != 0.0f);
  //Filter eventSelectionFilter = (evsel::sel8 == true);
  using filteredCollision = Filtered<Join<Collisions, EvSels>>::iterator;

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsMass{"nBinsMass", 100, "N bins in invariant mass histo"};
  Configurable<float> minPi0Mass{"minPi0Mass", 0.1f, "Maximum pi^0 Invariant Mass [GeV]"};
  Configurable<float> maxPi0Mass{"maxPi0Mass", 0.2f, "Maximum pi^0 Invariant Mass [GeV]"};
  Configurable<float> maxPi0Pt{"maxPi0Pt", 1.0f, "Maximum pi^0 Transverse Momentum [GeV]"};

  HistogramRegistry histosPi0{"histosPi0", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&) {
    const AxisSpec axisPi0Mass{nBinsMass, minPi0Mass, maxPi0Mass, "M_{#gamma#gamma_{PCM}} [GeV]"};
    const AxisSpec axisPi0Pt{nBinsPt, 0.0f, maxPi0Pt, "p_{T, candidate} [GeV]"};
    const AxisSpec axisAlpha{100, -1.0f, 1.0f, "#alpha"};
    const AxisSpec axisQt{100, 0.0f, 0.25f, "q_{T} [GeV]"};

    histosPi0.add("Pi0PCM", "#pi^{0} Candidates From Conversion Photons", kTH2F, {axisPi0Mass, axisPi0Pt});
    histosPi0.add("Pi0PCMArmenterosPodolanski", "Armenteros-Podolanski Plot for #pi^{0} Candidates (PCM)", kTH2F, {axisAlpha, axisQt});
  }

  void process(filteredCollision const&, ConversionPhoton const& photons1, LambdaHyperon const& photons2) {
    for (auto const& [photon1, photon2]: combinations(CombinationsStrictlyUpperIndexPolicy(photons1, photons2))) {
      TLorentzVector const photon1LorentzVector(photon1.px(), photon1.py(), photon1.pz(), photon1.e());
      TLorentzVector const photon2LorentzVector(photon2.px(), photon2.py(), photon2.pz(), photon2.e());
      TLorentzVector const pi0LorentzVector = photon1LorentzVector + photon2LorentzVector;
      float const pi0Mass = pi0LorentzVector.M();
      float const pi0Pt = pi0LorentzVector.Pt();

      TVector3 photon1Momentum(photon1.px(), photon1.py(), photon1.pz());
      TVector3 photon2Momentum(photon2.px(), photon2.py(), photon2.pz());
      TVector3 pi0Momentum(pi0LorentzVector.Px(), pi0LorentzVector.Py(), pi0LorentzVector.Pz());
      float pLongPos = photon1Momentum.Dot(pi0Momentum)/pi0Momentum.Mag();
      float pLongNeg = photon2Momentum.Dot(pi0Momentum)/pi0Momentum.Mag();
      float alpha = (pLongPos - pLongNeg)/(pLongPos + pLongNeg);
      float qT = photon1Momentum.Perp(pi0Momentum);

      if (pi0Mass < maxPi0Mass && pi0Mass > minPi0Mass) {
        histosPi0.fill(HIST("Pi0PCM"), pi0Mass, pi0Pt);
        histosPi0.fill(HIST("Pi0PCMArmenterosPodolanski"), alpha, qT);
      }
    }
  }
};

struct ReconstructMCSigma0viaPCM {
  Configurable<float> zVertexCut{"zVertexCut", 10.0f, "Maximum Primary Vertex Z coordinate [cm]"};
  Filter zVertexFilter = (nabs(collision::posZ) < zVertexCut);
  Filter zVertexErrorFilter = (collision::posZ != 0.0f);
  //Filter eventSelectionFilter = (evsel::sel8 == true);
  using filteredCollision = Filtered<Join<Collisions, EvSels>>::iterator;

  Configurable<int> nBinsPt{"nBinsPt", 120, "N bins in pT histo"};
  Configurable<int> nBinsMass{"nBinsMass", 150, "N bins in invariant mass histo"};
  Configurable<float> minSigma0Mass{"minSigma0Mass", 1.12f, "Maximum Sigma^0 Invariant Mass [GeV]"};
  Configurable<float> maxSigma0Mass{"maxSigma0Mass", 1.32f, "Maximum Sigma^0 Invariant Mass [GeV]"};
  Configurable<float> maxSigma0Pt{"maxSigma0Pt", 12.0f, "Maximum Sigma^0 Hyperon Transverse Momentum [GeV]"};

  HistogramRegistry histosMCSigma0{"histosMCSigma0", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&) {
    const AxisSpec axisSigma0Mass{nBinsMass, minSigma0Mass, maxSigma0Mass, "M_{#Lambda#gamma_{PCM}} [GeV]"};
    const AxisSpec axisSigma0Pt{nBinsPt, 0.0f, maxSigma0Pt, "p_{T, candidate} [GeV]"};

    histosMCSigma0.add("MCSigma0PCMgenLambdaGenPhoton", "#Sigma^{0} Candidates From Conversion Photons (gen. #Lambda, gen. #gamma)", kTH2F, {axisSigma0Mass, axisSigma0Pt});
    histosMCSigma0.add("MCSigma0PCMgenLambdaRecPhoton", "#Sigma^{0} Candidates From Conversion Photons (gen. #Lambda, rec. #gamma)", kTH2F, {axisSigma0Mass, axisSigma0Pt});
    histosMCSigma0.add("MCSigma0PCMrecLambdaGenPhoton", "#Sigma^{0} Candidates From Conversion Photons (rec. #Lambda, gen. #gamma)", kTH2F, {axisSigma0Mass, axisSigma0Pt});
    histosMCSigma0.add("MCSigma0PCMrecLambdaRecPhoton", "#Sigma^{0} Candidates From Conversion Photons (rec. #Lambda, rec. #gamma)", kTH2F, {axisSigma0Mass, axisSigma0Pt}); 
  }

  void process(filteredCollision const&, mcConversionPhoton const& mcPhotons, mcLambdaHyperon const& mcLambdas) {
    for (auto const& [mcLambda, mcPhoton]: combinations(CombinationsFullIndexPolicy(mcLambdas, mcPhotons))) {
      TLorentzVector const generatedLambdaLorentzVector(mcLambda.genpx(), mcLambda.genpy(), mcLambda.genpz(), mcLambda.gene());
      TLorentzVector const generatedPhotonLorentzVector(mcPhoton.genpx(), mcPhoton.genpy(), mcPhoton.genpz(), mcPhoton.gene());
     
      TLorentzVector const reconstructedLambdaLorentzVector(mcLambda.recpx(), mcLambda.recpy(), mcLambda.recpz(), mcLambda.rece());
      TLorentzVector const reconstructedPhotonLorentzVector(mcPhoton.recpx(), mcPhoton.recpy(), mcPhoton.recpz(), mcPhoton.rece());

      TLorentzVector const generatedLambdaGeneratedPhoton = generatedLambdaLorentzVector + generatedPhotonLorentzVector;
      TLorentzVector const generatedLambdaReconstructedPhoton = generatedLambdaLorentzVector + reconstructedPhotonLorentzVector;
      TLorentzVector const reconstructedLambdaGeneratedPhoton = reconstructedLambdaLorentzVector + generatedPhotonLorentzVector;
      TLorentzVector const reconstructedLambdaReconstructedPhoton = reconstructedLambdaLorentzVector + reconstructedPhotonLorentzVector;

      float Sigma0Mass = generatedLambdaGeneratedPhoton.M();
      float Sigma0Pt = generatedLambdaGeneratedPhoton.Pt();
      if (Sigma0Mass < maxSigma0Mass && Sigma0Mass > minSigma0Mass) {
        histosMCSigma0.fill(HIST("MCSigma0PCMgenLambdaGenPhoton"), Sigma0Mass, Sigma0Pt);           
      }

      Sigma0Mass = generatedLambdaReconstructedPhoton.M();
      Sigma0Pt = generatedLambdaReconstructedPhoton.Pt();
      if (Sigma0Mass < maxSigma0Mass && Sigma0Mass > minSigma0Mass) {
        histosMCSigma0.fill(HIST("MCSigma0PCMgenLambdaRecPhoton"), Sigma0Mass, Sigma0Pt);           
      }

      Sigma0Mass = reconstructedLambdaGeneratedPhoton.M();
      Sigma0Pt = reconstructedLambdaGeneratedPhoton.Pt();
      if (Sigma0Mass < maxSigma0Mass && Sigma0Mass > minSigma0Mass) {
        histosMCSigma0.fill(HIST("MCSigma0PCMrecLambdaGenPhoton"), Sigma0Mass, Sigma0Pt);           
      }

      Sigma0Mass = reconstructedLambdaReconstructedPhoton.M();
      Sigma0Pt = reconstructedLambdaReconstructedPhoton.Pt();
      if (Sigma0Mass < maxSigma0Mass && Sigma0Mass > minSigma0Mass) {
        histosMCSigma0.fill(HIST("MCSigma0PCMrecLambdaRecPhoton"), Sigma0Mass, Sigma0Pt);           
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
    adaptAnalysisTask<ReconstructSigma0viaPCM>(cfgc),
    adaptAnalysisTask<ReconstructNeutralPionsViaPCM>(cfgc),
    adaptAnalysisTask<ReconstructMCSigma0viaPCM>(cfgc)
  };
}