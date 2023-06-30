/*
 * RunAnalysis.C
 *
 *  Created on: Jun 17, 2019
 *      Author: mbenoit
 */

{

    gROOT->ProcessLine(".include ../../../src");
    gROOT->ProcessLine(".L ../../../lib/libAllpixObjects.so");
    gROOT->ProcessLine(".L ../../../lib/libAllpixCore.so");
    gROOT->ProcessLine(".L  ../../../lib/libAllpixModuleDetectorHistogrammer.so");
    gROOT->ProcessLine(".L Analysis.C++g");

    
    AnalysisB10("/home/esilva_sta/allpix-squared/examples/neutron_detection/analysis/",
		"output_ThermalNeutrons_Timepix3.root",
            "B10CoatedDetector");

    
    // AnalysisPE("/home/imb/Projects/allpix-squared/examples/neutron_detection/FastNeutronResults/",
    //         "output_FastNeutrons_measurement_Timepix3.root",
    //         "PECoatedDetector");

}
