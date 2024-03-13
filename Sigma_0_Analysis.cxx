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
/// \since 20/09/2023

#include "Framework/AnalysisTask.h" // needed in any case
#include "Framework/runDataProcessing.h" // needed in any case
#include "Framework/ASoA.h" // For columns and tables
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

struct Sigma0Reconstruction {
  float electronMass = o2::constants::physics::MassElectron;
  float chargedPionMass = o2::constants::physics::MassPionCharged;
  float protonMass = o2::constants::physics::MassProton;

  // General:
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsMass{"nBinsMass", 100, "N bins in invariant mass histo"};
  Configurable<float> etaCut{"etaCut", 1.2f, "Maximum Pseudorapidity"};

  // Collisions-related:
  Configurable<float> zVertexCut{"zVertexCut", 10.0f, "Maximum Primary Vertex Z coordinate [cm]"};
  Filter zVertexFilter = nabs(collision::posZ) < zVertexCut;
  Filter eventSelectionFilter = evsel::sel8 == true;
  using filteredCollision = Filtered<Join<Collisions, EvSels>>::iterator;

  // Decay vertex-related:
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1.0f, "DCA V0 Daughters [cm]"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06f, "DCA Pos To PV [cm]"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06f, "DCA Neg To PV [cm]"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.0f, "V0 Radius [cm]"};

  Configurable<float> maxPhotonAlpha{"maxPhotonAlpha", 0.9f, "Maximum Photon Decay Asymmetry"};
  Configurable<float> maxPhotonQt{"maxPhotonQt", 0.9f, "Maximum Photon Decay q_{T} [GeV]"};

  Configurable<float> maxLambdaAlpha{"maxLambdaAlpha", 0.99f, "Maximum Lambda Decay Asymmetry"};
  Configurable<float> minLambdaAlpha{"minLambdaAlpha", 0.0f, "Minimum Lambda Decay Asymmetry"};
  Configurable<float> maxLambdaQt{"maxLambdaQt", 0.9f, "Maximum Lambda Decay q_{T} [GeV]"};
  Configurable<float> minLambdaQt{"minLambdaQt", 0.0f, "Minimum Lambda Decay  [GeV]"};

  // filter can only be applied to static columns
  Filter preV0Filter = nabs(v0data::dcapostopv) > v0setting_dcapostopv && nabs(v0data::dcanegtopv) > v0setting_dcanegtopv && v0data::dcaV0daughters < v0setting_dcav0dau;
  using filteredV0s = Filtered<Join<V0Datas, McV0Labels>>;

  // Tracks-related:
  Configurable<float> maxPhotonMass{"maxPhotonMass", 0.1f, "Maximum Electron-Positron Invariant Mass [GeV]"};
  Configurable<float> maxPhotonPt{"maxPhotonPt", 0.5f, "Maximum Photon Transverse Momentum [GeV]"};

  Configurable<float> minLambdaMass{"minLambdaMass", 1.10f, "Minimum Lambda Hyperon Invariant Mass [GeV]"};
  Configurable<float> maxLambdaMass{"maxLambdaMass", 1.13f, "Maximum Lambda Hyperon Invariant Mass [GeV]"};
  Configurable<float> maxLambdaPt{"maxLambdaPt", 5.0f, "Maximum Lambda Hyperon Transverse Momentum [GeV]"};

  Configurable<float> maxTpcNSigmaEl{"maxTpcNSigmaEl", 3.0f, "Maximum N_{#sigma_{e}} from TPC signal"};
  Configurable<float> maxTpcNSigmaPr{"maxTpcNSigmaPr", 3.0f, "Maximum N_{#sigma_{p}} from TPC signal"};
  Configurable<float> maxTpcNSigmaPi{"maxTpcNSigmaPi", 3.0f, "Maximum N_{#sigma_{#pi}} from TPC signal"};

  Configurable<float> maxTrackDCA{"maxTrackDCA", 1.0f, "Maximum Track DCA [cm]"};

  Filter etaFilter = nabs(track::eta) < etaCut;
  Filter dcaFilter = nabs(track::dcaXY) < maxTrackDCA;

  using photonDaughterTracks = Join<Tracks, TracksDCA, pidTPCEl, McTrackLabels>;
  using filteredPhotonDaughterTracks = Filtered<photonDaughterTracks>;
  using lambdaDaughterTracks = Join<Tracks, TracksDCA, pidTPCPi, pidTPCPr, McTrackLabels>;
  using filteredLambdaDaughterTracks = Filtered<lambdaDaughterTracks>;

  HistogramRegistry histosCollisions{"histosCollisions", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosMC{"histosMC", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosConversionPhoton{"histosConversionPhoton", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosLambda{"histosLambda", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&) {
    // define axes:
    const AxisSpec axisVertexZ{100, -20.0f, 20.0f, "Vertex Z Coordinate [cm]"};
    const AxisSpec axisPhotonMass{nBinsMass, 0.0f, maxPhotonMass, "M_{e} [GeV]"};
    const AxisSpec axisLambdaMass{nBinsMass, minLambdaMass, maxLambdaMass, "M_{#Lambda} [GeV]"};
    const AxisSpec axisPhotonPt{nBinsPt, 0.0f, maxPhotonPt, "p_{T} [GeV]"};
    const AxisSpec axisLambdaPt{nBinsPt, 0.0f, maxLambdaPt, "p_{T} [GeV]"};
    const AxisSpec axisEta{150, -1.5f, 1.5f, "#eta"};
    const AxisSpec axisDCA{100, -5.0f, 5.0f, "DCA [cm]"};
    const AxisSpec axisV0Radius{100, 0.0f, 50.0f, "V0 Radius [cm]"};
    const AxisSpec axisAlpha{100, -1.0f, 1.0f, "#alpha"};
    const AxisSpec axisQt{100, 0.0f, 1.0f, "q_{T} [GeV]"};

    // collision-related histograms:
    histosCollisions.add("vertexZ", "vertexZ", kTH1F, {axisVertexZ});
    histosCollisions.add("eventCounter", "eventCounter", kTH1F, {{1, 0, 1, " Accepted Events Count"}});

    // Monte-Carlo histograms:
    histosMC.add("mcGenEta", "mcGenEta", kTH1F, {axisEta});
    histosMC.add("mcGenPt", "mcGenPt", kTH1F, {{100, 0, 5, "p_{T} [GeV]"}});
    histosMC.add("mcGenLambdaPt", "mcGenLambdaPt", kTH1F, {axisLambdaPt});
    histosMC.add("mcGenSigmaPt", "mcGenSigmaPt", kTH1F, {axisLambdaPt});
    histosMC.add("mcGenPhotonPt", "mcGenPhotonPt", kTH1F, {axisPhotonPt});
    histosMC.add("mcGenPhotonFromSigmaPt", "mcGenPhotonFromSigmaPt", kTH1F, {axisPhotonPt});
    histosMC.add("mcGenLambdaFromSigmaPt", "mcGenLambdaFromSigmaPt", kTH1F, {axisLambdaPt});

    // conversion photon-related histograms:
    histosConversionPhoton.add("photonEta", "photonEta", kTH1F, {axisEta});
    histosConversionPhoton.add("photonMass", "photonMass", kTH1F, {axisPhotonMass});
    histosConversionPhoton.add("photonPt", "photonPt", kTH1F, {axisPhotonPt});
    histosConversionPhoton.add("posTPCEl", "posTPCEl", kTH1F, {{100, -5, 5, "N_{#sigma_{e}}"}});
    histosConversionPhoton.add("negTPCEl", "negTPCEl", kTH1F, {{100, -5, 5, "N_{#sigma_{e}}"}});
    histosConversionPhoton.add("photonArmenterosPodolanski", "photonArmenterosPodolanski", kTH2F, {axisAlpha, axisQt}); 
    histosConversionPhoton.add("DCAv0Daughters", "DCAv0Daughters", kTH1F, {{100, 0.0f, 1.0f, "DCA [cm]"}});
    histosConversionPhoton.add("DCAv0Pos", "DCAv0Pos", kTH1F, {axisDCA});
    histosConversionPhoton.add("DCAv0Neg", "DCAv0Neg", kTH1F, {axisDCA});
    histosConversionPhoton.add("DCAv0cosPA", "DCAv0cosPA", kTH1F, {{100, 0.9f, 1.0f, "cos#theta"}});
    histosConversionPhoton.add("DCAv0Radius", "DCAv0Radius", kTH1F, {axisV0Radius});

    histosConversionPhoton.add("RecGenPhotonEnergyDifference", "RecGenPhotonEnergyDifference", kTH1F, {{100, -0.5f, 0.5f, "E_{rec} - E_{gen} [GeV]"}}); 
    histosConversionPhoton.add("RecGenPhotondEvsGenE", "RecGenPhotondEvsGenE", kTH2F, {{100, -0.5f, 0.5f, "E_{rec} - E_{gen} [GeV]"}, axisPhotonPt}); 
    histosConversionPhoton.add("MismatchPhoton", "MismatchPhoton", kTH1F, {{100, 0, maxPhotonPt, "p_{T} [GeV]"}});
    histosLambda.add("RecGenPhotonPtDifference", "RecGenPhotonPtDifference", kTH1F, {{100, -0.5, 0.5, "p_{T, rec} - p_{T, gen} [GeV]"}}); 
    histosLambda.add("RecGenPhotondPtvsGenPt", "RecGenPhotondPtvsGenPt", kTH2F, {{100, -0.5, 0.5, "p_{T, rec} - p_{T, gen} [GeV]"}, axisPhotonPt}); 

    // lambda-related histograms:
    histosLambda.add("LambdaEta", "LambdaEta", kTH1F, {axisEta});
    histosLambda.add("LambdaMass", "LambdaMass", kTH1F, {axisLambdaMass});
    histosLambda.add("LambdaPt", "LambdaPt", kTH1F, {axisLambdaPt});
    histosLambda.add("posTPCPr", "posTPCPr", kTH1F, {{100, -5, 5, "N_{#sigma_{p}}"}});
    histosLambda.add("negTPCPr", "negTPCPr", kTH1F, {{100, -5, 5, "N_{#sigma_{p}}"}});     
    histosLambda.add("posTPCPi", "posTPCPi", kTH1F, {{100, -5, 5, "N_{#sigma_{#pi}}"}});
    histosLambda.add("negTPCPi", "negTPCPi", kTH1F, {{100, -5, 5, "N_{#sigma_{#pi}}"}}); 
    histosLambda.add("LambdaArmenterosPodolanski", "LambdaArmenterosPodolanski", kTH2F, {axisAlpha, axisQt}); 
    histosLambda.add("DCAv0Daughters", "DCAv0Daughters", kTH1F, {{100, 0.0f, 1.0f, "DCA [cm]"}});
    histosLambda.add("DCAv0Pos", "DCAv0Pos", kTH1F, {axisDCA});
    histosLambda.add("DCAv0Neg", "DCAv0Neg", kTH1F, {axisDCA});
    histosLambda.add("DCAv0cosPA", "DCAv0cosPA", kTH1F, {{100, 0.9f, 1.0f, "cos#theta"}});
    histosLambda.add("DCAv0Radius", "DCAv0Radius", kTH1F, {axisV0Radius});

    histosLambda.add("RecGenLambdaEnergyDifference", "RecGenLambdaEnergyDifference", kTH1F, {{100, -0.5, 0.5, "E_{rec} - E_{gen} [GeV]"}}); 
    histosLambda.add("RecGenLambdadEvsGenE", "RecGenLambdadEvsGenE", kTH2F, {{100, -0.5, 0.5, "E_{rec} - E_{gen} [GeV]"}, axisLambdaPt}); 
    histosLambda.add("MismatchLambda", "MismatchLambda", kTH1F, {{100, 0, maxLambdaPt, "p_{T} [GeV]"}});
    histosLambda.add("RecGenLambdaPtDifference", "RecGenLambdaPtDifference", kTH1F, {{100, -0.5, 0.5, "p_{T, rec} - p_{T, gen} [GeV]"}}); 
    histosLambda.add("RecGenLambdadPtvsGenPt", "RecGenLambdadPtvsGenPt", kTH2F, {{100, -0.5, 0.5, "p_{T, rec} - p_{T, gen} [GeV]"}, axisLambdaPt}); 
  }

  void processCollisions(Collision const& collision) {
    histosCollisions.fill(HIST("vertexZ"), collision.posZ());
    if (abs(collision.posZ()) < zVertexCut) histosCollisions.fill(HIST("eventCounter"), 0.5);
  }
  PROCESS_SWITCH(Sigma0Reconstruction, processCollisions, "Process general collisions info", true);

  void processGeneratedCollisions(McCollision const& collision, McParticles const& mcparticles) {
    for (auto const& mcparticle: mcparticles) {
      histosMC.fill(HIST("mcGenEta"), mcparticle.eta());
      histosMC.fill(HIST("mcGenPt"), mcparticle.pt());
      if (abs(mcparticle.pdgCode()) == 3122) {
        histosMC.fill(HIST("mcGenLambdaPt"), mcparticle.pt());
        auto const& mother = mcparticle.mothers_first_as<McParticles>();
        if (abs(mother.pdgCode()) == 3212) histosMC.fill(HIST("mcGenLambdaFromSigmaPt"), mcparticle.pt());
      }

      if (mcparticle.pdgCode() == 22) {
        histosMC.fill(HIST("mcGenPhotonPt"), mcparticle.pt());
        auto const& mother = mcparticle.mothers_first_as<McParticles>();
        if (abs(mother.pdgCode()) == 3212) histosMC.fill(HIST("mcGenPhotonFromSigmaPt"), mcparticle.pt());
      }

      if (abs(mcparticle.pdgCode()) == 3212) {
        histosMC.fill(HIST("mcGenSigmaPt"), mcparticle.pt());
      }
    }
  }
  PROCESS_SWITCH(Sigma0Reconstruction, processGeneratedCollisions, "Process MC collisions info", true);

  void processConversionPhotons(filteredCollision const& collision, filteredV0s const& v0s, filteredPhotonDaughterTracks const&, McParticles const&) {
    for (auto const& v0: v0s) {
      // "Filter" on dynamic columns
      if (v0.v0cosPA() < v0setting_cospa) continue;
      if (v0.v0radius() < v0setting_radius) continue;
      if (abs(v0.alpha()) > maxPhotonAlpha) continue;
      if (v0.qtarm() > maxPhotonQt) continue;

      auto const& posPhotonDaughterTrack = v0.posTrack_as<filteredPhotonDaughterTracks>();
      auto const& negPhotonDaughterTrack = v0.negTrack_as<filteredPhotonDaughterTracks>();

      float const tpcNPosSigmaEl = posPhotonDaughterTrack.tpcNSigmaEl();
      float const tpcNNegSigmaEl = negPhotonDaughterTrack.tpcNSigmaEl();

      if (abs(tpcNPosSigmaEl) > maxTpcNSigmaEl || abs(tpcNNegSigmaEl) > maxTpcNSigmaEl) continue;

      float photonPx = posPhotonDaughterTrack.px() + negPhotonDaughterTrack.px();
      float photonPy = posPhotonDaughterTrack.py() + negPhotonDaughterTrack.py();
      float photonPz = posPhotonDaughterTrack.pz() + negPhotonDaughterTrack.pz();
      float electronEnergy = std::sqrt(negPhotonDaughterTrack.p()*negPhotonDaughterTrack.p() + electronMass*electronMass);
      float positronEnergy = std::sqrt(posPhotonDaughterTrack.p()*posPhotonDaughterTrack.p() + electronMass*electronMass);
      float photonEnergy = electronEnergy + positronEnergy;
      TLorentzVector photon(photonPx, photonPy, photonPz, photonEnergy);
      if (photon.Pt() > maxPhotonPt) continue;
      if (photon.M() > maxPhotonMass || photon.M() < 0) continue;

      histosConversionPhoton.fill(HIST("posTPCEl"), tpcNPosSigmaEl);
      histosConversionPhoton.fill(HIST("negTPCEl"), tpcNNegSigmaEl);

      histosConversionPhoton.fill(HIST("DCAv0Daughters"), v0.dcaV0daughters());
      histosConversionPhoton.fill(HIST("DCAv0Pos"), v0.dcapostopv());
      histosConversionPhoton.fill(HIST("DCAv0Neg"), v0.dcanegtopv());
      histosConversionPhoton.fill(HIST("DCAv0cosPA"), v0.v0cosPA());
      histosConversionPhoton.fill(HIST("DCAv0Radius"), v0.v0radius());
      histosConversionPhoton.fill(HIST("photonArmenterosPodolanski"), v0.alpha(), v0.qtarm());

      histosConversionPhoton.fill(HIST("photonMass"), photon.M());
      histosConversionPhoton.fill(HIST("photonPt"), photon.Pt());
      histosConversionPhoton.fill(HIST("photonEta"), photon.Eta());

      if (v0.has_mcParticle()) {
        auto const& v0mcParticle = v0.mcParticle();
        // check that the V0 comes from a photon
        if (v0mcParticle.pdgCode() == 22) {
          if (v0mcParticle.has_mothers()) {
            auto const& photonMother = v0mcParticle.mothers_first_as<McParticles>();
            // check that photon comes from Sigma^0 hyperon:
            if (abs(photonMother.pdgCode()) == 3212) {
              histosConversionPhoton.fill(HIST("RecGenPhotonEnergyDifference"), photonEnergy - v0mcParticle.e());
              histosConversionPhoton.fill(HIST("RecGenPhotondEvsGenE"), photonEnergy - v0mcParticle.e(), v0mcParticle.e());
              histosConversionPhoton.fill(HIST("RecGenPhotonPtDifference"), photon.Pt() - v0mcParticle.pt());
              histosConversionPhoton.fill(HIST("RecGenPhotondPtvsGenPt"), photon.Pt() - v0mcParticle.pt(), v0mcParticle.pt());
            } else {
              histosConversionPhoton.fill(HIST("MismatchPhoton"), v0mcParticle.pt());    
            }
          }
        }
      }
    }   
  }
  PROCESS_SWITCH(Sigma0Reconstruction, processConversionPhotons, "Process conversion photons", true);

  void processLambdas(filteredCollision const& collision, filteredV0s const& v0s, filteredLambdaDaughterTracks const&, McParticles const&) {
    for (auto const& v0: v0s) {
      // "Filter" on dynamic columns
      if (v0.v0cosPA() < v0setting_cospa) continue;
      if (v0.v0radius() < v0setting_radius) continue;
      if (abs(v0.alpha()) > maxLambdaAlpha || abs(v0.alpha()) < minLambdaAlpha) continue;
      if (v0.qtarm() > maxLambdaQt || v0.qtarm() < minLambdaQt) continue;

      auto const& posLambdaDaughterTrack = v0.posTrack_as<filteredLambdaDaughterTracks>();
      auto const& negLambdaDaughterTrack = v0.negTrack_as<filteredLambdaDaughterTracks>();

      float const tpcNPosSigmaPr = posLambdaDaughterTrack.tpcNSigmaPr();
      float const tpcNPosSigmaPi = posLambdaDaughterTrack.tpcNSigmaPi();

      float const tpcNNegSigmaPr = negLambdaDaughterTrack.tpcNSigmaPr();
      float const tpcNNegSigmaPi = negLambdaDaughterTrack.tpcNSigmaPi();

      bool posTraskIsProton = abs(tpcNPosSigmaPr) < maxTpcNSigmaPr;
      bool posTraskIsPion = abs(tpcNPosSigmaPr) < maxTpcNSigmaPi;

      bool negTraskIsProton = abs(tpcNNegSigmaPr) < maxTpcNSigmaPr;
      bool negTraskIsPion = abs(tpcNNegSigmaPi) < maxTpcNSigmaPi;

      bool isLambda = false;
      bool isAntiLambda = false;

      if (posTraskIsProton && negTraskIsPion) isLambda = true;
      if (negTraskIsProton && posTraskIsPion) isAntiLambda = true;

      if (!isLambda && !isAntiLambda) continue;

      float LambdaPx = posLambdaDaughterTrack.px() + negLambdaDaughterTrack.px();
      float LambdaPy = posLambdaDaughterTrack.py() + negLambdaDaughterTrack.py();
      float LambdaPz = posLambdaDaughterTrack.pz() + negLambdaDaughterTrack.pz();
      float protonEnergy = 0;
      float pionEnergy = 0;
      if (isLambda) {
        protonEnergy = std::sqrt(posLambdaDaughterTrack.p()*posLambdaDaughterTrack.p() + protonMass*protonMass);
        pionEnergy = std::sqrt(negLambdaDaughterTrack.p()*negLambdaDaughterTrack.p() + chargedPionMass*chargedPionMass);
      } else if (isAntiLambda) {
        protonEnergy = std::sqrt(negLambdaDaughterTrack.p()*negLambdaDaughterTrack.p() + protonMass*protonMass);
        pionEnergy = std::sqrt(posLambdaDaughterTrack.p()*posLambdaDaughterTrack.p() + chargedPionMass*chargedPionMass);
      }
      float LambdaEnergy = protonEnergy + pionEnergy;
      TLorentzVector Lambda(LambdaPx, LambdaPy, LambdaPz, LambdaEnergy);
      if (Lambda.Pt() > maxLambdaPt) continue;
      if (Lambda.M() > maxLambdaMass || Lambda.M() < minLambdaMass) continue;

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
      histosLambda.fill(HIST("DCAv0cosPA"), v0.v0cosPA());
      histosLambda.fill(HIST("DCAv0Radius"), v0.v0radius());
      histosLambda.fill(HIST("LambdaArmenterosPodolanski"), v0.alpha(), v0.qtarm());

      histosLambda.fill(HIST("LambdaMass"), Lambda.M());
      histosLambda.fill(HIST("LambdaPt"), Lambda.Pt());
      histosLambda.fill(HIST("LambdaEta"), Lambda.Eta());

      if (v0.has_mcParticle()) {
        auto const& v0mcParticle = v0.mcParticle();
        // check that the V0 comes from a lambda hyperon
        if (abs(v0mcParticle.pdgCode()) == 3122) {
          if (v0mcParticle.has_mothers()) {
            auto const& lambdaMother = v0mcParticle.mothers_first_as<McParticles>();
            // check that lambda hyperon comes from Sigma^0 hyperon:
            if (abs(lambdaMother.pdgCode()) == 3212) { 
              histosLambda.fill(HIST("RecGenLambdaEnergyDifference"), LambdaEnergy - v0mcParticle.e());
              histosLambda.fill(HIST("RecGenLambdadEvsGenE"), LambdaEnergy - v0mcParticle.e(), v0mcParticle.e());
              histosLambda.fill(HIST("RecGenLambdaPtDifference"), Lambda.Pt() - v0mcParticle.pt());
              histosLambda.fill(HIST("RecGenLambdadPtvsGenPt"), Lambda.Pt() - v0mcParticle.pt(), v0mcParticle.pt());
            } else {
              histosLambda.fill(HIST("MismatchLambda"), v0mcParticle.pt());       
            }
          }
        }
      }
    }   
  }
  PROCESS_SWITCH(Sigma0Reconstruction, processLambdas, "Process Lambda Hyperons", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{
    adaptAnalysisTask<Sigma0Reconstruction>(cfgc)
  };
}