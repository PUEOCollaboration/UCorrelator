#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all namespaces;

#pragma link C++ namespace pueo+;
#pragma link C++ namespace pueo::UCorrelator+;
#pragma link C++ namespace pueo::UCorrelator::gui+;
#pragma link C++ namespace pueo::UCorrelator::flags;
#pragma link C++ namespace pueo::UCorrelator::peakfinder;
#pragma link C++ namespace pueo::UCorrelator::shape;
#pragma link C++ namespace pueo::UCorrelator::spectrum;
#pragma link C++ namespace pueo::UCorrelator::image;

#pragma link C++ class pueo::UCorrelator::KDE2D;
#pragma link C++ class pueo::UCorrelator::peakfinder::FineMaximum;
#pragma link C++ class pueo::UCorrelator::peakfinder::RoughMaximum;
#pragma link C++ class pueo::UCorrelator::AntennaPositions;
#pragma link C++ class pueo::UCorrelator::Baseline;
#pragma link C++ class pueo::UCorrelator::TimeDependentAverage;
#pragma link C++ class pueo::UCorrelator::TimeDependentAverageLoader;
#pragma link C++ class pueo::UCorrelator::Correlator;
#pragma link C++ class pueo::UCorrelator::ComplicatedNotchFilter;
#pragma link C++ class pueo::UCorrelator::SineSubtractFilter;
#pragma link C++ class pueo::UCorrelator::AdaptiveFilterAbby;
#pragma link C++ class pueo::UCorrelator::CombinedSineSubtractFilter;
#pragma link C++ class pueo::UCorrelator::AdaptiveMinimumPhaseFilter;
#pragma link C++ class pueo::UCorrelator::AdaptiveBrickWallFilter; 
#pragma link C++ class pueo::UCorrelator::AdaptiveButterworthFilter;
#pragma link C++ class pueo::UCorrelator::WaveformCombiner;
#pragma link C++ class pueo::UCorrelator::Analyzer;
#pragma link C++ class pueo::UCorrelator::AnalysisConfig;
#pragma link C++ class pueo::UCorrelator::NoiseMachine;
#pragma link C++ class pueo::UCorrelator::PointingResolution;
#pragma link C++ class pueo::UCorrelator::ProbabilityMap;
#pragma link C++ class pueo::UCorrelator::ProbabilityMap::Params;
#pragma link C++ class pueo::UCorrelator::ProbabilityMap::Params::CollisionDetectionParams;
#pragma link C++ class pueo::UCorrelator::ProbabilityMap::Params::BackwardParams;
#pragma link C++ class pueo::UCorrelator::ProbabilityMap::Params::MCParams;
#pragma link C++ class pueo::UCorrelator::PointingResolutionModel;
#pragma link C++ class pueo::UCorrelator::PointingResolutionModelPlusHeadingError;
#pragma link C++ class pueo::UCorrelator::PointingResolutionParSNRModel;
#pragma link C++ class pueo::UCorrelator::ConstantPointingResolutionModel;
#pragma link C++ class pueo::UCorrelator::HeadingErrorEstimator;
#pragma link C++ class pueo::UCorrelator::EASFitter;
#pragma link C++ class pueo::UCorrelator::EASFitResult;
#pragma link C++ class pueo::UCorrelator::gui::Map+; 
#pragma link C++ class pueo::UCorrelator::gui::SummaryText; 



#endif

