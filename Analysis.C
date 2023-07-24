/*
 * Analysis.C
 *
 *  Created on: Jun 15, 2019
 *      Author: mbenoit
 */
#pragma includepath "../../../src:../../../src/objects";

#include <Math/DisplacementVector2D.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <stdio.h>
#include <fstream>
// FIXME: these includes should be absolute and provided with installation?
#include "../../../src/objects/MCParticle.hpp"
#include "../../../src/objects/Pixel.hpp"
#include "../../../src/objects/PixelCharge.hpp"
#include "../../../src/objects/PixelHit.hpp"
#include "../../../src/objects/PropagatedCharge.hpp"
#include "../../../src/modules/DetectorHistogrammer/Cluster.hpp"

#ifdef __MAKECINT__
#pragma link C++ class vector < allpix::PixelHit* > +;
#pragma link C++ class vector < allpix::MCTrack* > +;
#pragma link C++ class vector < allpix::MCParticle* > +;
#endif

int max_charge = 0;
int max_posX = -1;
int max_posY = -1;
int histo_iter = 0;

std::string format(const char* formatchar, ...) {
    va_list args;
    va_start(args, formatchar);
    std::string format(formatchar);
    size_t len = std::vsnprintf(NULL, 0, format.c_str(), args);
    va_end(args);
    std::vector<char> vec(len + 1);
    va_start(args, formatchar);
    std::vsnprintf(&vec[0], len + 1, format.c_str(), args);
    va_end(args);
    return &vec[0];
}

std::vector<allpix::Cluster> doClustering(std::vector<allpix::PixelHit*> hits) {
    std::vector<allpix::Cluster> clusters;
    std::map<const allpix::PixelHit*, bool> usedPixel;

    
    if(hits.size() < 1) {
        return clusters;
    }

    auto pixel_it = hits.begin();
    for(; pixel_it != hits.end(); pixel_it++) {
        const allpix::PixelHit* pixel_hit = (*pixel_it);

        // Check if the pixel has been used:
        if(usedPixel[pixel_hit]) {
            continue;
        }

        // Create new cluster
        allpix::Cluster cluster(pixel_hit);
        usedPixel[pixel_hit] = true;
	    //cout << "Creating new cluster with seed: " << pixel_hit->getPixel().getIndex().X()<< " " <<
        //pixel_hit->getPixel().getIndex().Y() <<std::endl;

        auto touching = [&](const allpix::PixelHit* pixel) {
         auto pxi1 = pixel->getIndex();
	        for(auto& cluster_pixel : cluster.getPixelHits()) {
                auto distance = [](unsigned int lhs, unsigned int rhs) { return (lhs > rhs ? lhs - rhs : rhs - lhs); };
                auto pxi2 = cluster_pixel->getIndex();
                if(distance(pxi1.x(), pxi2.x()) <= 1 && distance(pxi1.y(), pxi2.y()) <= 1) {
                    return true;
                }
            }
            return false;
        };

        // Keep adding pixels to the cluster:
        for(auto other_pixel = pixel_it + 1; other_pixel != hits.end(); other_pixel++) {
            const allpix::PixelHit* neighbor = (*other_pixel);

	        // Check if neighbor has been used or if it touches the current cluster:
            if(usedPixel[neighbor] || !touching(neighbor)) {
                continue;
            }
            
            
            cluster.addPixelHit(neighbor);
	        //cout  << "Adding pixel: " << neighbor->getPixel().getIndex().X() << " " << neighbor->getPixel().getIndex().Y()
	        //<< " Pixel Charge: " << neighbor->getSignal()
	        //<< std::endl;


	    
			if (neighbor->getSignal()>max_charge) {
			  max_charge = neighbor->getSignal();
			  max_posX = neighbor->getPixel().getIndex().X();
			  max_posY = neighbor->getPixel().getIndex().Y();
			}
		       
            usedPixel[neighbor] = true;
            other_pixel = pixel_it;
            

        }

	  //cout << " Max charge =  " << max_charge;
      //cout << " Min charge =  " << min_charge;
	  //cout << " Max posX : " << max_posX << " Max posY : " << max_posY << endl;
	  //return max_charge;
	
     clusters.push_back(cluster);

    }
    return clusters;
}

const allpix::PixelHit* FindSeedPixel(allpix::Cluster& cluster) {

    int temp_tot = 0;
    int max_index = 0;
    int idx = 0;
    const allpix::PixelHit* seed_pixel = new allpix::PixelHit();

    for(auto pixel : cluster.getPixelHits()) {
        if(pixel->getSignal() > temp_tot) {
            temp_tot = pixel->getSignal();
            max_index = idx;
            seed_pixel = pixel;
        }
        idx++;
     }

   
   //std::cout << format("%i, %i, ", seed_pixel->getIndex().x(), seed_pixel->getIndex().y());
    
    return seed_pixel;


}

void TrimCluster(allpix::Cluster& cluster, int N = 10) {

    auto seed_pixel = FindSeedPixel(cluster);
    auto new_cluster = allpix::Cluster(seed_pixel);
    auto seed_pixel_idx = seed_pixel->getIndex();

    for(auto pixel : cluster.getPixelHits()) {
        auto pix_addr = pixel->getIndex();
        if((pix_addr.x() == seed_pixel_idx.x()) && pix_addr.y() == seed_pixel_idx.y()) {

        } else if(abs(int((pix_addr.x() - seed_pixel_idx.x()))) <= N && abs(int(pix_addr.y() - seed_pixel_idx.y())) <= N) {
            new_cluster.addPixelHit(pixel);
        } else {
          //std::cout << format("Trimming pixel X: %i Y:%i", pix_addr.x(), pix_addr.y()) << std::endl;
        }
    }

    cluster = new_cluster;
}

/*
void PrintCluster(allpix::Cluster& cluster) {
      std::cout << "--------- CLUSTER begin -----------" << std::endl;

      for(auto pixel : cluster.getPixelHits()) {
        auto pix_addr = pixel->getIndex();
	//cluster_->SetBinContent(pix_addr.x(), pix_addr.y(), getSignal());
	//std::cout << format("X : %i Y : %i charge : %d",  pix_addr.x(), pix_addr.y(), pixel->getSignal()) << std::endl;

	  
	cout
	  //<< "X: " << pix_addr.x() << "Y: " << pix_addr.y()
	  << "  " << pixel->getSignal();
       
        ;
		
	//loop through pixels in the cluster


      }	
    std::cout << "--------- CLUSTER end   -----------" << std::endl;
    
}
*/  

/*
void PrintClusterMap(allpix::Cluster& cluster) {
  std::cout << "--------- CLUSTER MAP begins -----------" << std::endl;
  for(auto pixel : cluster.getPixelHits()) {
    auto pix_addr = pixel->getIndex();
    // histo_iter++;
    // if (histo_iter==1) cluster_1->Fill(pix_addr.x(), pix_addr.y(), pixel->getSignal());
      //else if (histo_iter==2) cluster_2->Fill(pix_addr.x(), pix_addr.y(), 1.0, pixel->getSignal());
      //else if (histo_iter==3) cluster_3->Fill(pix_addr.x(), pix_addr.y(), 1.0, pixel->getSignal()); 
    //else if (histo_iter==4) cluster_4->SetBinContent(p.ix_addr.x(), pix_addr.y(), getSignal());
    //else if (histo_iter==5) cluster_5->SetBinContent(pix_addr.x(), pix_addr.y(), getSignal());
    //else if (histo_iter==6) cluster_6->SetBinContent(pix_addr.x(), pix_addr.y(), getSignal());
    // cluster_->SetBinContent(pix_addr.x(), pix_addr.y(), getSignal());
    std::cout << format("X : %i Y : %i ", pix_addr.x(), pix_addr.y()) << std::endl;
    ;
  }

  std::cout << "--------- CLUSTER MAP ends   -----------" << std::endl;
}
*/

/* 
void PrintPHit(allpix::PixelHit& pixelHit) {

  std::cout << "--------- pixelHit begin -----------" << std::endl;
  for(auto pixel : pixelHit.getPixelHits()) {
    auto pix_addr = pixel->getIndex();
    std::cout << format("X : %i Y : %i ", pix_addr.x(), pix_addr.y()) << std::endl;
    ;
  }
  std::cout << "--------- pixelHit end   -----------" << std::endl;
}
*/

void AnalysisB10(string input_file_folder, 
            string input_file_name,
            string detector_name,
            double pitch = 0.055,
            int mat_size = 256,
            int maxTOT = 2e9,
            int chargeTOT = 400e6) {

    //gROOT->ProcessLine(".L /eda/allpix2/allpix-squared/lib/libAllpixObjects.so");
    //gROOT->ProcessLine(".L /eda/allpix2/allpix-squared/lib/libAllpixModuleDetectorHistogrammer.so");
    
    //cout << "Opening file " << input_file_name << std::endl;
    auto file = TFile::Open((input_file_folder + "/" +input_file_name).c_str());
    //cout << "file opened" << std::endl;

    TTree* pixel_hit_tree = static_cast<TTree*>(file->Get("PixelHit"));
    TTree* mctrack_tree = static_cast<TTree*>(file->Get("MCTrack"));
    TTree* mcparticles_tree = static_cast<TTree*>(file->Get("MCParticle")); 
    int nevents = mcparticles_tree->GetEntries();
    //cout << "got the trees" << std::endl;

    //cout << "branches picked" << std::endl;

    std::vector<allpix::PixelHit*> input_hits_bot;
    std::vector<allpix::MCTrack*> input_tracks;
    std::vector<allpix::MCParticle*> input_particles_bot;
    //cout << "created vectors" << std::endl;

    pixel_hit_tree->FindBranch(detector_name.c_str())->SetObject(&input_hits_bot);
    mctrack_tree->FindBranch("global")->SetObject(&input_tracks);
    mcparticles_tree->FindBranch(detector_name.c_str())->SetObject(&input_particles_bot);
    //cout << "Assigned object vectors to branches" << std::endl;

    std::string output_file = input_file_folder +"/"+ detector_name + "_output_plots.root";
    TFile* outfile = new TFile(output_file.c_str(), "recreate");
    outfile->cd();

    //Booking Histograms
    TH1D* resx_bot = new TH1D("resx_bot", "Residual X  ", 200, -0.05, 0.05);
    TH1D* resy_bot = new TH1D("resy_bot", "Residual Y ", 200, -0.05, 0.05);

    TH1D* resx_point = new TH1D("resx_point", "Local Points Difference in X", 200, -0.01, 0.01);
    TH1D* resy_point = new TH1D("resy_point", "Local Points Difference in Y", 200, -0.01, 0.01);

    TH1D* resx_max = new TH1D("resx_max", "Max Charge Residual X  ", 200, -0.05, 0.05);
    TH1D* resy_max = new TH1D("resy_max", "Max Charge Residual Y  ", 200, -0.05, 0.05);
    
    TH1I* clu_tot_bot = new TH1I("clu_tot_bot", "Cluster TOT", 1000, 0, maxTOT);    

    TH2D* clusize_charge = new TH2D("clusize_charge", "Cluster Charge vs Cluster Size", 100, 0, 100, 100, 0, maxTOT);

    TH1D* hitmapx_bot = new TH1D("hitmapx_bot", "Cluster position in X", 2500, 0, 2*7.04);

    TH1D* hitmapy_bot = new TH1D("hitmapy_bot", "Cluster position in Y", 2500, 0, 2*7.04);

    TH1I* clu_size = new TH1I("clu_size", "Cluster Size", 50, 0, 50);

    TH2D* hitmap_bot = new TH2D("hitmap_bot", "Cluster position", 500, 0, 2*7.04, 500, 0, 2*7.04);

    TH1D* pix_maxcharge = new TH1D("pix_maxcharge", "Max Charge", 100, 0, chargeTOT);

    TH1D* maxX_position = new TH1D("maxX_position", "Max Charge Position in X", 257, -1, 255 );

    TH1D* maxY_position = new TH1D("maxY_position", "Max Charge Position in Y", 257, -1, 255);    
        
    TH2D* cluster_1 = new TH2D("cluster_1", "Cluster1 Pixels and Charge", 256, 0, 257, 256, 0, 257);
    TH2D* cluster_2 = new TH2D("cluster_2", "Cluster2 Pixels and Charge", 256, 0, 257, 256, 0, 257);
  
    TH1D* hitmapx_mc = new TH1D("hitmapx_mc", "MC Position in X", 2500, 0, 2*7.04);
    TH1D* hitmapy_mc = new TH1D("hitmapy_mc", "MC Position in Y", 2500, 0, 2*7.04);

    //histogram of track distance sqrt[(McX - McX_end)^2 + (McY - McY_end)^2]
    TH1D* distance = new TH1D("distance", "Track Distance", 1000, 0, 0.007);

    //histogram of the Z axis end point
    TH1D* z_endpt = new TH1D("z_endpt", "Z-axis End Point", 1000, -0.151, -0.143);

    //histogram of alpha angles (degrees)
    TH1D* est_alpha_angle = new TH1D("est_alpha_angle", "Estimated Cluster Alpha Angles", 1000, 0, 100);
    //histogram of beta angles (degrees)
    TH1D* beta_angle = new TH1D("beta_angle", "Estimated Beta Angles", 1000, -100, 100);


    //variables for next portion of code

    int n_alpha = 0;
    int n_lithium = 0;
    int n_conversion = 0;
    int n_coinc = 0;
    int n_SiCapture = 0;

    
    double resx_bot_coinc = 0;
    double resy_bot_coinc = 0;

    // csv header
    //std::cout << " Particle type, cluster charge, cluster position in x, cluster position in y, local Start X, local start Y, local start Z, local end X, local end Y, local end Z, X, Y, charge " << std::endl; 
    
    for(int i = 0; i < pixel_hit_tree->GetEntries(); ++i) {
        int totalTOT = 0;
        //quit after some number of events for testing:
        //if (i > 1000) continue; 
	    //cout << "----------------------Event " << i << "---------------------" << std::endl;
	    //cout << "################### MCTruth #######################" << std::endl;
        pixel_hit_tree->GetEvent(i);
        mctrack_tree->GetEvent(i);
        mcparticles_tree->GetEvent(i);
        bool has_int_bot = false;
        bool has_int_top = false;
        for(auto& p : input_tracks) {
            int id = p->getParticleID();
            auto start = p->getStartPoint();
            auto end = p->getEndPoint();
            auto energy = p->getKineticEnergyInitial() - p->getKineticEnergyFinal();
	    
            //cout << format("PID: %i ", id) << endl;
	        //cout << format("(%f %f %f ) (%f %f %f ) dedx : %f MeV",
	             //                           start.X(),
			     //                           start.Y(),
			     //                           start.Z(),
			     //                           end.X(),
			     //                           end.Y(),
			     //                           end.Z(),
			     //		                      energy)
	             //                           << std::endl;
            if(id == 1000140290 || id == 1000140300) {
                n_SiCapture++;
            }
        }
        //cout << "################### Bottom sensor ################" << std::endl;
	
     max_charge = 0;
     max_posX = -1;
     max_posY = -1;
     
     auto clusters = doClustering(input_hits_bot);
	
     //if (clusters.size()>0) cout << format("Found %i clusters, with %i pixels", clusters.size(), input_hits_bot.size()) << std::endl;

	        for(auto& p : input_particles_bot) {
                //cout << format("PID: %i ", p->getParticleID()) << endl;
		    }

		
	 auto comp_clu = [](allpix::Cluster& a, allpix::Cluster& b) { return a.getCharge() > b.getCharge(); };
        std::sort(clusters.begin(), clusters.end(), comp_clu);
        for(auto& clu : clusters) {
             //std::cout << "old cluster" << std::endl;
	         //if (clusters.size == 2) {
	         // PrintCluster(clu);
			 // PrintClusterMap(clu);
			 //}

               	TrimCluster(clu);
            
			 //std::cout << "new cluster" << std::endl;
			 //PrintCluster(clu);

             double posX = clu.getPosition().X();
             double posY = clu.getPosition().Y();
	      

	     int cluTOT = clu.getCharge();
	     //if (cluTOT>0) cout << format(" X : %f Y : %f  %i, size %i ", clu.getPosition().X(), clu.getPosition().Y(), cluTOT, clu.getSize()) << std::endl;
	    
            for(auto p : clu.getMCParticles()) {

                if(p == NULL) {
                    continue;
                }
                int id = p->getParticleID();
                
                if(((id == 1000030070) || (id == 1000020040)) && has_int_bot == false){
		        //if((id == 1000030070) && has_int_bot == false) {
    	        //if((id == 1000020040) && has_int_bot == false) { 
                    double mcX_g = p->getGlobalStartPoint().X();
                    double mcY_g = p->getGlobalStartPoint().Y();
                    double mcX = p->getLocalStartPoint().X();
                    double mcY = p->getLocalStartPoint().Y(); 
                    double mcZ = p->getLocalStartPoint().Z();                  
		            //double mcX_g_end = p->getGlobalEndPoint().X();
		            //double mcY_g_end = p->getGlobalEndPoint().Y();
		            double mcX_end = p->getLocalEndPoint().X();
		            double mcY_end = p->getLocalEndPoint().Y();
                    double mcZ_end = p->getLocalEndPoint().Z();

		            double pi = M_PI;
                    double delta_y = mcY_end - mcY;
                    double delta_x = mcX_end - mcX;
                    double distances =  sqrt(((mcY_end - mcY) * (mcY_end - mcY)) + ((mcX_end - mcX) * (mcX_end - mcX)));
                    double z_len = abs(mcZ_end - mcZ);
                    double angle_deg = atan(distances/z_len)*(180/pi);
                    double beta_angle_deg = atan(delta_y/delta_x)*(180/pi);
                    
                    

                    cout  << " alpha angle is " << angle_deg << " & beta angle is " << beta_angle_deg << ::endl;
		    
		         // cout << "Neutron interaction detected" << std::endl;

                    if(cluTOT > 0 && clu.getSize() > 1) {

			         resx_point->Fill(mcX_end - mcX);
			         resy_point->Fill(mcY_end - mcY);

			         resx_max->Fill((max_posX*0.055) - mcX);
			         resy_max->Fill((max_posY*0.055) - mcY);
			
			         resx_bot->Fill(posX - mcX);
                     resy_bot->Fill(posY - mcY);

			         clu_size->Fill(clu.getSize());
			         clusize_charge->Fill(clu.getSize(), cluTOT);
			         pix_maxcharge->Fill(max_charge);
			         maxX_position->Fill(max_posX);
			         maxY_position->Fill(max_posY);
            
			         distance->Fill(distances);
                     z_endpt->Fill(mcZ_end);
                     est_alpha_angle->Fill(angle_deg);
                     beta_angle->Fill(beta_angle_deg);
                   
                     // const allpix::PixelHit* pixel_hit = (*pixel_it);

                     /*
                     for (i = 0, i < 256, i++) {
                     for (j = 0, j < 256, j++)
                     pixel_hit.i.x()
                     //set pixel to
                     auto pxi_addr = pixel_hit->getIndex();
                     cout << " " << pixel_hit->getSignal();
                     }
                     */


			         //cout << "posX : " << posX << " mc truth X : " << mcX << endl;
			         //cout << "posY : " << posY << " mc truth Y : " << mcY << endl;
            
                     //int i = 0;
			         //double clu_array[256]; 

                       
			            for(auto pixel : clu.getPixelHits()) {
			             auto pix_addr = pixel->getIndex();
			             // histo_iter++;
			             if (histo_iter==1) cluster_1->Fill(pix_addr.x(), pix_addr.y(), pixel->getSignal()); 
			             else if (histo_iter==2) cluster_2->Fill(pix_addr.x(), pix_addr.y(), pixel->getSignal());
			              //else if (histo_iter==3) cluster_3->Fill(pix_addr.x(), pix_addr.y(), 1.0, pixel->getSignal());
			             //else if (histo_iter==4) cluster_4->SetBinContent(pix_addr.x(), pix_addr.y(), getSignal());
			             //else if (histo_iter==5) cluster_5->SetBinContent(pix_addr.x(), pix_addr.y(), getSignal());
			             //else if (histo_iter==6) cluster_6->SetBinContent(pix_addr.x(), pix_addr.y(), getSignal());
			             // cluster_->SetBinContent(pix_addr.x(), pix_addr.y(), getSignal());
			             // std::cout << format("X : %i Y : %i ", pix_addr.x(), pix_addr.y()) << std::endl;
 			             ;
	
    			        }
			         /*  
                     histo_iter++;

                     int clu_event_iter = 0;
                     bool printed = false;
                     //std::cout << " i = " << i << std::endl;
      
                        for(auto pixel : clu.getPixelHits()) {
			                auto pix_addr = pixel->getIndex();
			                
                            //if (clu_event_iter < 10) 
                            if (clu_event_iter == 0 && !printed) {
                                if (id == 1000030070) std::cout << "0, ";
                                else if(id == 1000020040) std::cout << "1, ";
                                std::cout << cluTOT << ", " << posX << ", " << posY << ", " << mcX << ", " << mcY << ", " << mcZ << ", " << mcX_end << ", " << mcY_end << ", " << mcZ_end;
                                printed = true; 
                                //break;seed_pixel->getIndex().x(), seed_pixel->getIndex().y()
                            }

			                if (clu_event_iter == 0) {
                             if (printed)
                                std::cout << ", ";
                                std:: cout << pix_addr.x() << ", " << pix_addr.y() << ", "; 
                                std::streamsize prevPrecision = std::cout.precision();
                                std::cout << std::fixed << std::setprecision(0) << pixel->getSignal();
                                std::cout.precision(prevPrecision);
                            }
                            ;
                        }
                     cout << std::endl;
    	             clu_event_iter++;
                     //std::cout << clu_event_iter << std::endl;
			            */ 
                    }
                  resx_bot_coinc = posX - mcX + pitch / 2;
                  resy_bot_coinc = posY - mcY + pitch / 2;
		   
		         if(clu.getSize() > 3)
		            clu_tot_bot->Fill(cluTOT);
                    totalTOT += cluTOT;
                    hitmap_bot->Fill(posX, posY);
                    hitmapx_bot->Fill(posX);
		            hitmapy_bot->Fill(posY);
		            hitmapx_mc->Fill(mcX);
		            hitmapy_mc->Fill(mcY);
		     
                    n_conversion++;
                    has_int_bot = true;
                }
                if(id == 1000030070) {
                    n_lithium++;
		            //cout << "Lithium ion cnt " << n_lithium << std::endl;
                }
                if(id == 1000020040) {
                    n_alpha++;
                    //cout << "Alpha cnt  " << n_alpha << std::endl;
                }

                if(has_int_bot) {
                    break;
                }
            }
        }
    }

    /*
    double efficiency = 100 * double(n_conversion) / nevents;
        cout << format("[Interaction report] %i Lithium, %i Alpha, %i coincidence, %i conversion detected, %i Si Capture, %f %% "
                   "efficiency ",
                   n_lithium,
                   n_alpha,
                   n_coinc,
                   n_conversion,
                   n_SiCapture,
                   efficiency)
         << std::endl;
    */

    std::vector<TH1D*> resplots{resx_bot, resy_bot};
    for(auto& plot : resplots) {
        plot->SetLineWidth(2);
        plot->SetLineColor(kBlue);
        // plot->SetFillColor(kRed);
        // plot->SetFillStyle(3353);
        plot->GetXaxis()->SetTitle("X_{neutron} - X_{reco} [mm]");
        gStyle->SetOptFit(0011);
    }   

    
    TCanvas* can = new TCanvas();
    can->Draw();
    can->SetWindowSize(1400, 800);
    can->Divide(2, 1);
    can->cd(1);
    resx_bot->Draw();
    can->cd(2);
    resy_bot->Draw();

    can->Print((input_file_folder + "ResidualXY_Bottom.png").c_str());

    TCanvas* can2 = new TCanvas();
    can2->Draw();
    can2->SetWindowSize(1400, 800);
    can2->Divide(2, 1);
    can2->cd(1);
    resx_point->Draw();
    can2->cd(2);
    resy_point->Draw();
    resx_point->GetXaxis()->SetTitle("[mm]");
    resx_point->GetYaxis()->SetTitle("number of pixels");
    resy_point->GetXaxis()->SetTitle("[mm]");
    resy_point->GetYaxis()->SetTitle("number of pixels");

    can2->Update();
    can2->Print((input_file_folder + "ResidualXY_Local_Point.png").c_str());
  
    TCanvas* can3 = new TCanvas();
    can3->Draw();
    can3->SetWindowSize(1400, 800);
    // can3->Divide(2);
    // can3->cd(1);
    clu_tot_bot->SetFillColor(kRed);
    clu_tot_bot->SetLineColor(kRed);
    clu_tot_bot->SetFillStyle(3353);
    clu_tot_bot->SetLineWidth(2);
    clu_tot_bot->GetXaxis()->SetRangeUser(0, maxTOT);
    clu_tot_bot->GetXaxis()->SetTitle("TOT [A.U.]");
    clu_tot_bot->GetYaxis()->SetTitle("number of clusters");  
    clu_tot_bot->SetStats(0);
    clu_tot_bot->Draw();
    // can3->cd(2);
    // clu_tot_top->Draw();

    can3->Print((input_file_folder + "clusterTOT_Bottom.png").c_str());
    
    TCanvas* can4 = new TCanvas();
    can4->Draw();
    can4->SetWindowSize(1400, 800);
    clusize_charge->Draw("colz"); //colz for 2D histograms
    clusize_charge->GetXaxis()->SetTitle("size [pixels]");
    clusize_charge->GetYaxis()->SetTitle("charge [eV]");
    clusize_charge->GetZaxis()->SetTitle("number of clusters");

    can4->Update();
    can4->Print((input_file_folder + "clusize_charge.png").c_str());

    TCanvas* can5 = new TCanvas();
    can5->Draw();
    can5->SetWindowSize(1400, 800);
    hitmapx_bot->Draw();
    hitmapx_bot->GetXaxis()->SetTitle("[mm]");
    hitmapx_bot->GetYaxis()->SetTitle("number of clusters");

    can5->Update();
    can5->Print((input_file_folder + "hitmap_X_bottom.png").c_str());

    TCanvas* can6 = new TCanvas();
    can6->Draw();
    can6->SetWindowSize(1400, 800);
    hitmapy_bot->Draw("colz");
    hitmapy_bot->GetXaxis()->SetTitle("[mm]");
    hitmapy_bot->GetYaxis()->SetTitle("number of clusters");

    can6->Update();
    can6->Print((input_file_folder + "hitmap_Y_bot.png").c_str());

    TCanvas* can7 = new TCanvas();
    can7->Draw();
    can7->SetWindowSize(1400, 800);
    clu_size->Draw(); 
    clu_size->GetXaxis()->SetTitle("cluster size [pixels]");
    clu_size->GetYaxis()->SetTitle("number of clusters");

    can7->Update();    
    can7->Print((input_file_folder + "Cluster_size.png").c_str());

    TCanvas* can8 = new TCanvas();
    can8->Draw();
    can8->SetWindowSize(1400, 800);
    hitmap_bot->Draw();
    hitmap_bot->GetXaxis()->SetTitle("[mm]");
    hitmap_bot->GetYaxis()->SetTitle("[mm]");

    can8->Update();  
    can8->Print((input_file_folder + "hitmap_bot.png").c_str());

    TCanvas* can9 = new TCanvas();
    can9->Draw();
    can9->SetWindowSize(1400, 800);
    pix_maxcharge->Draw();
    pix_maxcharge->GetXaxis()->SetTitle("pixel charge [eV]");
    pix_maxcharge->GetYaxis()->SetTitle("number of clusters");

    can9->Update();  
    can9->Print((input_file_folder + "Max_Charge.png").c_str());

    TCanvas* can10 = new TCanvas();
    can10->Draw();
    can10->SetWindowSize(1400, 800);
    maxX_position->Draw();
    maxX_position->GetXaxis()->SetTitle("position [pixel]");
    maxX_position->GetYaxis()->SetTitle("number of clusters");

    can10->Update(); 
    can10->Print((input_file_folder + "Max_X_position.png").c_str());

    TCanvas* can11 = new TCanvas();
    can11->Draw();
    can11->SetWindowSize(1400, 800);
    can11->Divide(2, 1);
    can11->cd(1);
    resx_max->Draw();
    can11->cd(2);
    resy_max->Draw();
    resx_max->GetXaxis()->SetTitle("[mm]");
    resy_max->GetYaxis()->SetTitle("number of pixels");

    can11->Update(); 
    can11->Print((input_file_folder + "ResidualXY_Max.png").c_str());

    TCanvas* can12 = new TCanvas();
    can12->Draw();
    can12->SetWindowSize(1400, 800);
    hitmapx_mc->Draw();
    hitmapx_mc->GetXaxis()->SetTitle("[mm]");
    hitmapx_mc->GetYaxis()->SetTitle("number of clusters");

    can12->Update();
    can12->Print((input_file_folder + "hitmap_X_MC.png").c_str());

    TCanvas* can13 = new TCanvas();
    can13->Draw();
    can13->SetWindowSize(1400, 800);
    hitmapy_mc->Draw();
    hitmapy_mc->GetXaxis()->SetTitle("[mm]");
    hitmapy_mc->GetYaxis()->SetTitle("number of clusters");

    can13->Update();
    can13->Print((input_file_folder + "hitmap_Y_MC.png").c_str());

    TCanvas* can140 = new TCanvas();
    can140->Draw();
    can140->SetWindowSize(1400, 800);
    maxY_position->Draw();
    maxY_position->GetXaxis()->SetTitle("position [pixel]");
    maxY_position->GetYaxis()->SetTitle("number of clusters");

    can140->Update(); 
    can140->Print((input_file_folder + "Max_Y_position.png").c_str());
    
    TCanvas* can14 = new TCanvas();
    can14->Draw();
    can14->SetWindowSize(1400, 800);
    cluster_1->Draw("colz");
    cluster_1->GetXaxis()->SetTitle("position [pixel]");
    cluster_1->GetYaxis()->SetTitle("position [pixel]");
    cluster_1->GetZaxis()->SetTitle("charge [eV]");

    can14->Update();
    can14->Print((input_file_folder + "cluster_1.png").c_str());

    TCanvas* can15 = new TCanvas();
    can15->Draw();
    can15->SetWindowSize(1400, 800);
    cluster_2->Draw("colz");
    cluster_2->GetXaxis()->SetTitle("position [pixel]");
    cluster_2->GetYaxis()->SetTitle("position [pixel]");
    cluster_2->GetZaxis()->SetTitle("charge [eV]");

    can15->Update();
    can15->Print((input_file_folder + "cluster_2.png").c_str());

    
    TCanvas* can16 = new TCanvas();
    can16->Draw();
    can16->SetWindowSize(1400, 800);
    distance->Draw();
    distance->GetXaxis()->SetTitle("[mm]");
    distance->GetYaxis()->SetTitle("number of tracks");

    can16->Update();
    can16->Print((input_file_folder + "distance.png").c_str());
    
   TCanvas* can17 = new TCanvas();
    can17->Draw();
    can17->SetWindowSize(1400, 800);
    z_endpt->Draw();
    z_endpt->GetXaxis()->SetTitle("[mm]");
    z_endpt->GetYaxis()->SetTitle("number of z-points");

    can17->Update();
    can17->Print((input_file_folder + "z_endpoint.png").c_str());

 TCanvas* can18 = new TCanvas();
    can18->Draw();
    can18->SetWindowSize(1400, 800);
    est_alpha_angle->Draw();
    est_alpha_angle->GetXaxis()->SetTitle("[degrees]");
    est_alpha_angle->GetYaxis()->SetTitle("number of angles");

    can18->Update();
    can18->Print((input_file_folder + "est_alpha_angle.png").c_str());

 TCanvas* can19 = new TCanvas();
    can19->Draw();
    can19->SetWindowSize(1400, 800);
    beta_angle->Draw();
    beta_angle->GetXaxis()->SetTitle("[degrees]");
    beta_angle->GetYaxis()->SetTitle("number of angles");

    can19->Update();
    can19->Print((input_file_folder + "beta_angle.png").c_str());


    //Write out the histogram
    resx_bot->Write();
    resy_bot->Write();

    resx_point->Write();
    resy_point->Write();

    clu_tot_bot->Write();

    clusize_charge->Write();    

    hitmapx_bot->Write();

    hitmapy_bot->Write();

    hitmapx_mc->Write();

    hitmapy_mc->Write();
    
    clu_size->Write();

    hitmap_bot->Write();

    pix_maxcharge->Write();

    maxX_position->Write();

    resx_max->Write();
    resy_max->Write();

    maxY_position->Write();
    
    cluster_1->Write();
    cluster_2->Write();
    
    distance->Write();

    z_endpt->Write();
   
   est_alpha_angle->Write();
   beta_angle->Write();
  

    /*
    cluster_3->Write();
  
    cluster_4->Write();
    cluster_5->Write();
    cluster_6->Write();
    */
    
    //resx_bot_g->Write();
    //resy_bot_g->Write();
    
    //resx_top->Write();
    //resy_top->Write();
    //resx_com->Write();
    //resy_com->Write();
    
    //hitmap_top->Write()
    //hitmapx_top->Write();
    //clu_tot_top->Write();
    //clu_tot_com->Write();

    // outfile->Close();
    // file->Close();
      
    //  delete mcparticles_tree;
    //  delete mctrack_tree;
    //  delete file;
    //  delete outfile;
    //  delete pixel_hit_tree;
    //  delete resx_bot;
    //	delete resy_bot;
    //	delete resx_top;
    //	delete resy_top;
    //	delete resx_com;
    //	delete resy_com;
    //	delete clu_tot_top;
    //	delete clu_tot_bot;
    //	delete hitmap_bot;
    //	delete hitmap_top;
}


// AnalysisPE is not being run

void AnalysisPE(string input_file_folder,
              string input_file_name,
              string detector_name,
              double pitch = 0.055,
              int mat_size = 256,
		int maxTOT = 2e9) {

    // gROOT->ProcessLine(".L /eda/allpix2/allpix-squared/lib/libAllpixObjects.so");
    // gROOT->ProcessLine(".L /eda/allpix2/allpix-squared/lib/libAllpixModuleDetectorHistogrammer.so");

    // std::string input_file_folder = format("examples/UCN_Detection/test/");
    cout << "Opening file " << input_file_name << std::endl;
    auto file = TFile::Open((input_file_folder + "/" +input_file_name).c_str());
    cout << "file opened" << std::endl;

    TTree* pixel_hit_tree = static_cast<TTree*>(file->Get("PixelHit"));
    TTree* mctrack_tree = static_cast<TTree*>(file->Get("MCTrack"));
    TTree* mcparticles_tree = static_cast<TTree*>(file->Get("MCParticle"));
    int nevents = mcparticles_tree->GetEntries();
    cout << "got the trees" << std::endl;

    cout << "branches picked" << std::endl;

    std::vector<allpix::PixelHit*> input_hits_bot;
    std::vector<allpix::MCTrack*> input_tracks;
    std::vector<allpix::MCParticle*> input_particles_bot;
    cout << "created vectors" << std::endl;

    pixel_hit_tree->FindBranch(detector_name.c_str())->SetObject(&input_hits_bot);
    mctrack_tree->FindBranch("global")->SetObject(&input_tracks);
    mcparticles_tree->FindBranch(detector_name.c_str())->SetObject(&input_particles_bot);
    cout << "Assigned object vectors to branches" << std::endl;

    std::string output_file = input_file_folder +"/"+ detector_name + "_output_plots.root";
    TFile* outfile = new TFile(output_file.c_str(), "recreate");
    outfile->cd();

    TH1D* resx_bot = new TH1D("resx_bot", "Residual X ", 200, -0.05, 0.05);
    TH1D* resy_bot = new TH1D("resy_bot", "Residual Y ", 200, -0.05, 0.05);

    TH2D* hitmap_bot = new TH2D("hitmap_bot", "Cluster position", 500, 0, 2*7.04, 500, 0, 2*7.04);

    TH1D* hitmapx_bot = new TH1D("hitmapx_bot", "Cluster position in X", 2500, 0, 2*7.04);
    
    TH1I* clu_tot_bot = new TH1I("clu_tot_bot", "Cluster TOT", 1000, 0, maxTOT);

    TH1I* clu_size = new TH1I ("clu_size", "Cluster Size", 30, 0, 30); 
    
    int n_alpha = 0;
    int n_lithium = 0;
    int n_conversion = 0;
    int n_coinc = 0;
    int n_SiCapture = 0;

    double resx_bot_coinc = 0;
    double resy_bot_coinc = 0;

    for(int i = 0; i < pixel_hit_tree->GetEntries(); ++i) {
        int totalTOT = 0;
	cout << "----------------------Event " << i << "---------------------" << std::endl;
	cout << "################### MCTruth #######################" << std::endl;
	pixel_hit_tree->GetEvent(i);
	mctrack_tree->GetEvent(i);
	mcparticles_tree->GetEvent(i);
	bool has_int_bot = false;
	bool has_int_top = false;
	for(auto& p : input_tracks) {
	  int id = p->getParticleID();
	  auto start = p->getStartPoint();
	  auto end = p->getEndPoint();
	  auto energy = p->getKineticEnergyInitial() - p->getKineticEnergyFinal();
	  cout << format("PID: %i ", id) << endl;
	  cout << format("(%f %f %f ) (%f %f %f ) dedx : %f MeV",
			 start.X(),
			 start.Y(),
			 start.Z(),
			 end.X(),
			 end.Y(),
			 end.Z(),
			 energy)

                 << std::endl;
            if(id == 1000140290 || id == 1000140300) {
                n_SiCapture++;
            }
        }
        cout << "################### Bottom sensor ################" << std::endl;
	//	int max_charge = 0;
        auto clusters = doClustering(input_hits_bot);
        cout << format("Found %i clusters, with %i pixels", clusters.size(), input_hits_bot.size()) << std::endl;

        for(auto& p : input_particles_bot) {
            cout << format("PID: %i ", p->getParticleID()) << endl;
        }

        auto comp_clu = [](allpix::Cluster& a, allpix::Cluster& b) { return a.getCharge() > b.getCharge(); };
        std::sort(clusters.begin(), clusters.end(), comp_clu);
        for(auto& clu : clusters) {

                   	std::cout << "old cluster" << std::endl;
			//                   	PrintCluster(clu);
                   	TrimCluster(clu);
            
                   	std::cout << "new cluster" << std::endl;
			//	PrintCluster(clu);

            double posX = clu.getPosition().X() ;
            double posY = clu.getPosition().Y();

	      
            int cluTOT = clu.getCharge();
            cout << format(" X : %f Y : %f  %i, size %i", clu.getPosition().X(), clu.getPosition().Y(), cluTOT, clu.getSize()) << std::endl;

            for(auto p : clu.getMCParticles()) {

                if(p == NULL) {
                    continue;
                }
                int id = p->getParticleID();

                if(id == 2212) {
                    // double mcX = p->getGlobalStartPoint().X();
                    // double mcY = p->getGlobalStartPoint().Y();
                    double mcX = p->getLocalStartPoint().X();
                    double mcY = p->getLocalStartPoint().Y();                   
                    cout << "Neutron interaction detected" << std::endl;
                    if(cluTOT > 0 && clu.getSize() > 1) {
                        resx_bot->Fill(posX - mcX );
                        resy_bot->Fill(posY - mcY );
		  	clu_size->Fill(clu.getSize()); 
                        cout << "posX : " << posX << " mc truth X : " << mcX << endl;
                        cout << "posY : " << posY << " mc truth Y : " << mcY << endl;
			
                    }
                    resx_bot_coinc = posX - mcX + pitch / 2;
                    resy_bot_coinc = posY - mcY + pitch / 2;

                    if(clu.getSize() > 3)
                        clu_tot_bot->Fill(cluTOT);
                    totalTOT += cluTOT;
                    hitmap_bot->Fill(posX, posY);
                    hitmapx_bot->Fill(posX);
                    n_conversion++;
                    has_int_bot = true;
                }
                if(id == 2212) {
                    n_lithium++;
                    cout << "Lithium ion cnt " << n_lithium << std::endl;
                }
                if(id == 1000020040) {
                    n_alpha++;
                    cout << "Alpha cnt  " << n_alpha << std::endl;
                }

                if(has_int_bot) {
                    break;
                }
            }
        }
    }



    double efficiency = 100 * double(n_conversion) / nevents;
    
    cout << format("[Interaction report] %i Protons, %i Alpha, %i coincidence, %i conversion detected, %i Si Capture, %f %% "
                   "efficiency ",
                   n_lithium,
                   n_alpha,
                   n_coinc,
                   n_conversion,
                   n_SiCapture,
                   efficiency)
         << std::endl;

    std::vector<TH1D*> resplots{resx_bot, resy_bot};
    for(auto& plot : resplots) {
        plot->SetLineWidth(2);
        plot->SetLineColor(kBlue);
        // plot->SetFillColor(kRed);
        // plot->SetFillStyle(3353);
        plot->GetXaxis()->SetTitle("X_{neutron} - X_{reco} [mm]");
        gStyle->SetOptFit(0011);
    }

    TCanvas* can = new TCanvas();
    can->Draw();
    can->SetWindowSize(1400, 800);
    can->Divide(2, 1);
    can->cd(1);
    resx_bot->Draw();
    can->cd(2);
    resy_bot->Draw();

    TCanvas* can7 = new TCanvas();
    can7->Draw();
    can7->SetWindowSize(1400, 800);
    clu_size->Draw();

    //    can->cd(3);
    //    resx_top->Draw();
    //    can->cd(4);
    //    resy_top->Draw();
    // auto fitx = resx_bot->Fit("gaus", "S", "", -0.01, 0.01);
    // auto fity = resy_bot->Fit("gaus", "S", "", -0.01, 0.01);

    // cout << format("resolution X : %f um Y: %f um", fitx->Parameter(2) * 1000, fity->Parameter(2) * 1000) << std::endl;

    can->Print((input_file_folder + "ResidualXY_Bottom.png").c_str());


    TCanvas* can2 = new TCanvas();
    can2->Draw();
    //can2->Divide(2);
    //can2->cd(1);
    hitmap_bot->Draw("colz");
    //can2->cd(2);
   // hitmap_top->Draw("colz");

    TCanvas* can3 = new TCanvas();
    can3->Draw();
    can3->SetWindowSize(1400, 800);

    // can3->Divide(2);
    // can3->cd(1);
    clu_tot_bot->SetFillColor(kRed);
    clu_tot_bot->SetLineColor(kRed);
    clu_tot_bot->SetFillStyle(3353);
    clu_tot_bot->SetLineWidth(2);
    clu_tot_bot->GetXaxis()->SetRangeUser(0, maxTOT);
    clu_tot_bot->GetXaxis()->SetTitle("TOT [A. U.]");
    clu_tot_bot->SetStats(0);
    clu_tot_bot->Draw();
    // can3->cd(2);
    // clu_tot_top->Draw();
    can3->Print((input_file_folder + "clusterTOT_Bottom.png").c_str());

    // TCanvas* can6 = new TCanvas();
    // can6->Draw();
    // can6->SetWindowSize(1400, 800);

    // // can3->Divide(2);
    // // can3->cd(1);
    // clu_tot_top->SetFillColor(kRed);
    // clu_tot_top->SetLineColor(kRed);
    // clu_tot_top->SetFillStyle(3353);
    // clu_tot_top->SetLineWidth(2);
    // clu_tot_top->GetXaxis()->SetRangeUser(0, maxTOT);
    // clu_tot_top->GetXaxis()->SetTitle("TOT [A. U.]");
    // clu_tot_top->SetStats(0);
    // clu_tot_top->Draw();
    // // can3->cd(2);
    // // clu_tot_top->Draw();
    // can6->Print((input_file_folder + "clusterTOT_Top.png").c_str());

    // TCanvas* can7 = new TCanvas();
    // can7->Draw();
    // can7->SetWindowSize(1400, 800);

    // // can3->Divide(2);
    // // can3->cd(1);
    // clu_tot_com->SetFillColor(kRed);
    // clu_tot_com->SetLineColor(kRed);
    // clu_tot_com->SetFillStyle(3353);
    // clu_tot_com->SetLineWidth(2);
    // clu_tot_com->GetXaxis()->SetRangeUser(0, 2 * maxTOT);
    // clu_tot_com->GetXaxis()->SetTitle("TOT [A. U.]");
    // clu_tot_com->SetStats(0);
    // clu_tot_com->Draw();
    // // can3->cd(2);
    // // clu_tot_top->Draw();
    // can7->Print((input_file_folder + "clusterTOT_Combined.png").c_str());

    TCanvas* can8 = new TCanvas();
    can8->Draw();
    can8->SetWindowSize(1400, 800);
    hitmapx_bot->Draw();
    can8->Print((input_file_folder + "hitmap_X_bottom.png").c_str());
    
    // TCanvas* can9 = new TCanvas();
    // can9->Draw();
    // can9->SetWindowSize(1400, 800);
    // hitmapx_top->Draw();
    // can9->Print((input_file_folder + "hitmap_X_top.png").c_str());

    resx_bot->Write();
    resy_bot->Write();
    clu_size->Write();
    //resx_top->Write();
    //resy_top->Write();
    //resx_com->Write();
    //resy_com->Write();
    hitmap_bot->Write();
    //hitmap_top->Write();
    hitmapx_bot->Write();
   // hitmapx_top->Write();
    //clu_tot_top->Write();
    clu_tot_bot->Write();
    //clu_tot_com->Write();

    // outfile->Close();
    //    file->Close();

    //    delete mcparticles_tree;
    //    delete mctrack_tree;
    //    delete file;
    //    delete outfile;
    //    delete pixel_hit_tree;
    //    delete resx_bot;
    //	delete resy_bot;
    //	delete resx_top;
    //	delete resy_top;
    //	delete resx_com;
    //	delete resy_com;
    //	delete clu_tot_top;
    //	delete clu_tot_bot;
    //	delete hitmap_bot;
    //	delete hitmap_top;
}
