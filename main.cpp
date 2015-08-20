#include <ostream>
#include<istream>
#include<fstream>
#include <iostream>
#include <string>
//#include "matrix.h"
//#include "homoest.h"

#include <stdlib.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/flann/flann.hpp>
#include <opencv2/legacy/legacy.hpp>
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/nonfree/nonfree.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "svd.h"
#include "inter.h"

using namespace std;
using namespace cv;
#define CLIP2(minv, maxv, value) (min(maxv, max(minv, value))) 
void dibMar(Mat img_matches,Mat img_object,Mat img_scene,Mat res_obj,vector<Point2f> scene_corners){
    
  circle(img_scene,scene_corners[0],4,Scalar(255,255,255),CV_FILLED,8);  
  circle(img_scene,scene_corners[1],4,Scalar(255,255,255),CV_FILLED,8);  
  circle(img_scene,scene_corners[2],4,Scalar(255,255,255),CV_FILLED,8); 
  circle(img_scene,scene_corners[3],4,Scalar(255,255,255),CV_FILLED,8); 
  /*circle(img_matches,scene_corners[1],4,Scalar(255,255,255),CV_FILLED,8);  
  circle(img_matches,scene_corners[2],4,Scalar(255,255,255),CV_FILLED,8);  
  circle(img_matches,scene_corners[3],4,Scalar(255,255,255),CV_FILLED,8); */ 
    

    
 /* line( img_matches, scene_corners[0] + Point2f( img_object.cols, 0), scene_corners[1] + Point2f( img_object.cols, 0), Scalar(0, 255, 0), 4 );
  line( img_matches, scene_corners[1] + Point2f( img_object.cols, 0), scene_corners[2] + Point2f( img_object.cols, 0), Scalar( 0, 255, 0), 4 );
  line( img_matches, scene_corners[2] + Point2f( img_object.cols, 0), scene_corners[3] + Point2f( img_object.cols, 0), Scalar( 0, 255, 0), 4 );
  line( img_matches, scene_corners[3] + Point2f( img_object.cols, 0), scene_corners[0] + Point2f( img_object.cols, 0), Scalar( 0, 255, 0), 4 );
  */
  imshow( "coordenadas", img_scene );
  waitKey(0);

}
Point2f mult(Mat H,Point2f obj){
   Point2f res; 
    res.x=H.at<double>(0,0)*obj.x+H.at<double>(0,1)*obj.y+H.at<double>(0,2);
    res.y=H.at<double>(1,0)*obj.x+H.at<double>(1,1)*obj.y+H.at<double>(1,2);   
    double val=H.at<double>(2,0)*obj.x+H.at<double>(2,1)*obj.y+H.at<double>(2,2);   
    res.x=res.x/val;
    res.y=res.y/val;
  return res;
}
double  calcularError(Mat H,vector<Point2f> obj,vector<Point2f> scene){
  double suma=0;
  Point2f scene2;
  for(int i=0;i<obj.size();i++){
      scene2=mult(H,scene[i]);
      suma+=sqrt(pow(obj[i].x-scene2.x,2)+pow(obj[i].y-scene2.y,2));
  }
  return suma;
}
Mat homografia_0(Mat img_object,Mat img_matches,vector<Point2f> obj,vector<Point2f> scene){
  Mat H = findHomography( obj, scene,0);// CV_LMEDS // CV_RANSAC // el parametro cero para usar todos los puntos
  /*vector<Point2f> obj_corners(4);
  obj_corners[0] = cvPoint(0,0); 
  obj_corners[1] = cvPoint( img_object.cols, 0 );
  obj_corners[2] = cvPoint( img_object.cols, img_object.rows ); 
  obj_corners[3] = cvPoint( 0, img_object.rows );
  vector<Point2f> scene_corners(4);
  perspectiveTransform( obj_corners, scene_corners, H);
  dibMar(img_matches,img_object,scene_corners);
  */
  
  return H;
}
Mat homografia_RANSAC(Mat img_object,Mat img_matches,vector<Point2f> obj,vector<Point2f> scene){
  Mat H = findHomography( obj, scene,CV_RANSAC);// CV_LMEDS // CV_RANSAC // el parametro cero para usar todos los puntos
  /*vector<Point2f> obj_corners(4);
  obj_corners[0] = cvPoint(0,0); obj_corners[1] = cvPoint( img_object.cols, 0 );
  obj_corners[2] = cvPoint( img_object.cols, img_object.rows ); obj_corners[3] = cvPoint( 0, img_object.rows );
  vector<Point2f> scene_corners(4);
  perspectiveTransform( obj_corners, scene_corners, H);
  dibMar(img_matches,img_object,scene_corners);*/
  return H;  
}
Mat homografia_LMEDS(Mat img_object,Mat img_matches,vector<Point2f> obj,vector<Point2f> scene){
  Mat H = findHomography( obj, scene,CV_LMEDS);// CV_LMEDS // CV_RANSAC // el parametro cero para usar todos los puntos
  /*vector<Point2f> obj_corners(4);
  obj_corners[0] = cvPoint(0,0); obj_corners[1] = cvPoint( img_object.cols, 0 );
  obj_corners[2] = cvPoint( img_object.cols, img_object.rows ); obj_corners[3] = cvPoint( 0, img_object.rows );
  vector<Point2f> scene_corners(4);
  perspectiveTransform( obj_corners, scene_corners, H);
  dibMar(img_matches,img_object,scene_corners);*/
  return H;
}
void persTrasn(Mat img_object,Mat img_matches,Mat H){
  vector<Point2f> obj_corners(4);
  obj_corners[0] = cvPoint(0,0); 
  obj_corners[1] = cvPoint( img_object.cols, 0 );
  obj_corners[2] = cvPoint( img_object.cols, img_object.rows ); 
  obj_corners[3] = cvPoint( 0, img_object.rows );
  vector<Point2f> scene_corners(4);
  perspectiveTransform( obj_corners, scene_corners, H);
  //dibMar(img_matches,img_object,scene_corners);
}



Mat detecDesc(Mat img_object,Mat img_scene,vector<Point2f> &obj,vector<Point2f> &scene){

     SiftFeatureDetector detector(2000);

     std::vector<KeyPoint> keypoints_object, keypoints_scene;

     detector.detect(img_object, keypoints_object);
     detector.detect(img_scene, keypoints_scene);

     SiftDescriptorExtractor extractor;

     Mat descriptors_object, descriptors_scene;

     extractor.compute( img_object, keypoints_object, descriptors_object );
     extractor.compute( img_scene, keypoints_scene, descriptors_scene );

     BruteForceMatcher<L2<float>> matcher;
     std::vector< DMatch > matches;
     matcher.match( descriptors_object, descriptors_scene, matches );

     double max_dist = 0; double min_dist = 1000;

     //-- Quick calculation of max and min distances between keypoints
     for( int i = 0; i < descriptors_object.rows; i++ ){ 
       double dist = matches[i].distance;
       if( dist < min_dist ) min_dist = dist;
       if( dist > max_dist ) max_dist = dist;
     }

     //cout<<"-- Max dist :  \n"<< max_dist ;
    // cout<<"-- Min dist :  \n"<< min_dist;

     std::vector< DMatch > good_matches;
     
     for( int i = 0; i < descriptors_object.rows; i++ ){ 
		 if( matches[i].distance < 3*min_dist ){ 
			 good_matches.push_back( matches[i]); 
                         //cout<<matches[i].queryIdx<<", ";
		 }
     }
     cout<<"numero de matches"<<good_matches.size();

     Mat img_matches;
     drawMatches( img_object, keypoints_object, img_scene, keypoints_scene,
                  good_matches, img_matches, Scalar::all(-1), Scalar::all(-1),//-1 es = a random color de cir y linea
                  vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );

     for( int i = 0; i < good_matches.size(); i++ ){
       obj.push_back( keypoints_object[ good_matches[i].queryIdx ].pt );
      // cout<<keypoints_object[ good_matches[i].queryIdx ].pt<<", ";
       scene.push_back( keypoints_scene[ good_matches[i].trainIdx ].pt );
     }
     cout<<endl<<" numero de puntos de correspondencia :"<<good_matches.size()<<endl;
    std::string folder = "correspSIFTa";
	std::ofstream results("/home/system/NetBeansProjects/kanatani/correspSIFTa/corresp1.txt");
        cout<<"----------------------------------------------------"<<endl;
     results<<obj.size()<<endl;
     for(int i=0;i<good_matches.size();i++)
        results<<obj[i].y<<endl;     
     
     for(int i=0;i<good_matches.size();i++)
        results<<obj[i].x<<endl;
     
     for(int i=0;i<good_matches.size();i++)
        results<<scene[i].y<<endl;
     for(int i=0;i<good_matches.size();i++)
        results<<scene[i].x<<endl;
     results.close();
     cout<<"-------------------------------------------------------"<<endl;
     return img_matches;
}
Mat inversa(Mat A){
    Mat result=Mat(3,3, CV_64F, double(0));
double determinant =    +A.at<double>(0,0)*(A.at<double>(1,1)*A.at<double>(2,2)-A.at<double>(2,1)*A.at<double>(1,2))
                        -A.at<double>(0,1)*(A.at<double>(1,0)*A.at<double>(2,2)-A.at<double>(1,2)*A.at<double>(2,0))
                        +A.at<double>(0,2)*(A.at<double>(1,0)*A.at<double>(2,1)-A.at<double>(1,1)*A.at<double>(2,0));
    double invdet = 1/determinant;
    result.at<double>(0,0) =  (A.at<double>(1,1)*A.at<double>(2,2)-A.at<double>(2,1)*A.at<double>(1,2))*invdet;
    result.at<double>(1,0) = -(A.at<double>(0,1)*A.at<double>(2,2)-A.at<double>(0,2)*A.at<double>(2,1))*invdet;
    result.at<double>(2,0) =  (A.at<double>(0,1)*A.at<double>(1,2)-A.at<double>(0,2)*A.at<double>(1,1))*invdet;
    result.at<double>(0,1) = -(A.at<double>(1,0)*A.at<double>(2,2)-A.at<double>(1,2)*A.at<double>(2,0))*invdet;
    result.at<double>(1,1) =  (A.at<double>(0,0)*A.at<double>(2,2)-A.at<double>(0,2)*A.at<double>(2,0))*invdet;
    result.at<double>(2,1) = -(A.at<double>(0,0)*A.at<double>(1,2)-A.at<double>(1,0)*A.at<double>(0,2))*invdet;
    result.at<double>(0,2) =  (A.at<double>(1,0)*A.at<double>(2,1)-A.at<double>(2,0)*A.at<double>(1,1))*invdet;
    result.at<double>(1,2) = -(A.at<double>(0,0)*A.at<double>(2,1)-A.at<double>(2,0)*A.at<double>(0,1))*invdet;
    result.at<double>(2,2) =  (A.at<double>(0,0)*A.at<double>(1,1)-A.at<double>(1,0)*A.at<double>(0,1))*invdet;
    return result;
}
Mat multMat(Mat A,Mat B){
    Mat R=Mat(3,3, CV_64F, double(0));
    int x,y,i;
   for(y=0; y<3; y++)
    for(x=0; x<3; x++)
      for(i=0; i<3; i++)         // which == second->rows
	R.at<double>(x,y) += A.at<double>(x,i) * B.at<double>(i,y);
    return R;
}
Mat dlt2d(Mat img_object,Mat img_matches,vector<Point2f> obj,vector<Point2f> scene){
    Mat H = Mat(3,3, CV_64F, double(0));
    int i=0;
    int n=obj.size();
     CvMat *A = cvCreateMat(2*n, 9, CV_64FC1); 
     CvMat *U = cvCreateMat(2*n, 2*n, CV_64FC1); 
     CvMat *D = cvCreateMat(2*n, 9, CV_64FC1); 
     CvMat *V = cvCreateMat(9, 9, CV_64FC1); 
     cvZero(A); 
    for(i=0; i<n; i++){ 
        // 2*i row 
        cvmSet(A,2*i,3,-obj[i].x); 
        cvmSet(A,2*i,4,-obj[i].y); 
        cvmSet(A,2*i,5,-1); 
        cvmSet(A,2*i,6,scene[i].y*obj[i].x); 
        cvmSet(A,2*i,7,scene[i].y*obj[i].y); 
        cvmSet(A,2*i,8,scene[i].y); 
        // 2*i+1 row 
        cvmSet(A,2*i+1,0,obj[i].x); 
        cvmSet(A,2*i+1,1,obj[i].y); 
        cvmSet(A,2*i+1,2,1); 
        cvmSet(A,2*i+1,6,-scene[i].x*obj[i].x); 
        cvmSet(A,2*i+1,7,-scene[i].x*obj[i].y); 
        cvmSet(A,2*i+1,8,-scene[i].x); 
     } 
    // SVD 
     // The flags cause U and V to be returned transposed 
     // Therefore, in OpenCV, A = U^T D V 
     cvSVD(A, D, U, V, CV_SVD_U_T|CV_SVD_V_T); 
       // SVD(A, D, U, V, CV_SVD_U_T|CV_SVD_V_T); 
        for(int i=0;i<9;i++){
            H.at<double>(i/3,i%3)=cvmGet(V, 8, i);
        }

       // H.at<double>(0,2) = (-1)*abs(v[2][s]);
       // H.at<double>(1,2) = (-1)*abs(v[5][s]);        
   // cout<<"]La matriz H:"<<H;
        /*vector<Point2f> obj_corners(4);
        obj_corners[0] = cvPoint(0,0); obj_corners[1] = cvPoint( img_object.cols, 0 );
        obj_corners[2] = cvPoint( img_object.cols, img_object.rows ); obj_corners[3] = cvPoint( 0, img_object.rows );
        vector<Point2f> scene_corners(4);
        perspectiveTransform( obj_corners, scene_corners, H);
        dibMar(img_matches,img_object,scene_corners);*/
    return H;
}  
Mat dlt2DNorm(Mat img_object,Mat img_matches,vector<Point2f> obj,vector<Point2f> scene){
	Mat H = Mat(3,3, CV_64F, double(0));
        float **A, **v;
        float *w;
        int m,n;
        int i,j;
        m = 2*obj.size();
        n = 9;
        A= (float**)malloc(sizeof(float*)*m);
        v= (float**)malloc(sizeof(float*)*n);
        w= (float*)malloc(sizeof(float)*n);
        for(j=0;j<m;j++) 
          A[j]=(float*)malloc(sizeof(float)*n);
        
        for(j=0;j<n;j++) 
          v[j]=(float*)malloc(sizeof(float)*n);        
        
        double c1x,c1y,c2x,c2y,s1,s2;
        Mat t1,t2;
        t1=Mat(3,3, CV_64F, double(0));
        t2=Mat(3,3, CV_64F, double(0));
        double sum11=0.0,sum12=0.0,sum21=0.0,sum22=0.0;
        
	for(int k=0;k<obj.size();k++){
            sum11=sum11+double(obj[k].x);
            sum12=sum12+double(obj[k].y);
            sum21=sum21+double(scene[k].x);
            sum22=sum22+double(scene[k].y);
        }
        
        c1x=sum11/(double)n;
        c1y=sum12/(double)n;
        c2x=sum21/(double)n;
        c2y=sum22/(double)n;
        double dist1=0;
        double dist2=0;

        for(i=0;i<n;i++) {
           dist1 +=  (obj[i].x - c1x)*(obj[i].x - c1x) + (obj[i].y - c1y)*(obj[i].y - c1y);
           dist2 +=  (scene[i].x - c2x)*(scene[i].x- c2x) + (scene[i].y - c2y)*(scene[i].y - c2y);
        }
        
        s1 = sqrt(2.0*n/dist1);
        s2 = sqrt(2.0*n/dist2);
        //cout<<" s1:"<<s1<<", s2:"<<s2<<endl;
        t1.at<double>(0,0)=s1;//s1      0       0
        t1.at<double>(1,0)=0;//0        s1      0
        t1.at<double>(2,0)=0;//-clxs1 -clys1    1
        t1.at<double>(0,1)=0;
        t1.at<double>(1,1)=s1;
        t1.at<double>(2,1)=0;
        t1.at<double>(0,2)=-c1x*s1;
        t1.at<double>(1,2)=-c1y*s1;
        t1.at<double>(2,2)=1;
        t2.at<double>(0,0)=s2;
        t2.at<double>(1,0)=0;
        t2.at<double>(2,0)=0;
        t2.at<double>(0,1)=0;
        t2.at<double>(1,1)=s2;
        t2.at<double>(2,1)=0;
        t2.at<double>(0,2)=-c2x*s2;
        t2.at<double>(1,2)=-c2y*s2;
        t2.at<double>(2,2)=1;
        
        for(int k=0;k<n;k++){
            obj[k].x=s1*(obj[k].x-c1x);
            obj[k].y=s1*(obj[k].y-c1y);
            scene[k].x=s2*(scene[k].x-c2x);
            scene[k].y=s2*(scene[k].y-c2y);       
        }

	 for(i=0; i<2*n; i+=2){
       A[i][0] = 0; 
       A[i][1] = 0; 
       A[i][2] = 0;
       A[i][3] = (-1)* obj[i].x ;
       A[i][4] = (-1)*obj[i].y ; 
       A[i][5] = -1;
       A[i][6] = scene[i].y * obj[i].x;
       A[i][7] = scene[i].y *obj[i].y;
       A[i][8] = scene[i].y;


       A[(i+1)][0] = obj[i].x;
       A[(i+1)][1] = obj[i].y; 
       A[(i+1)][2] = 1;
       A[(i+1)][3] = 0; 
       A[(i+1)][4] = 0; 
       A[(i+1)][5] = 0;
       A[(i+1)][6] = (-1)*scene[i].x * obj[i].x;
       A[(i+1)][7] = (-1)*scene[i].x* obj[i].y;
       A[(i+1)][8] = (-1)*scene[i].x; 
    }
 
       int r=dsvd(A, 2*n, 9, w, v);
        double minsingularvalue = 9999999;
        int s = -1;        
        for(i=0;i<n;i++) {           
          if (w[i] < minsingularvalue) {
            minsingularvalue = w[i];
            s = i;
          }
          //cout<<w[i]<<endl;
        } 
         /*for(int i=0;i<9;i++){
            H.at<double>(i/3,i%3)=abs(v[i][s]);
        }*/
       H.at<double>(0,0) = abs(v[0][s]);
        H.at<double>(1,0) = abs(v[3][s]);
        H.at<double>(2,0) = abs(v[6][s]);
        H.at<double>(0,1) = abs(v[1][s]);
        H.at<double>(1,1) = abs(v[4][s]);
        H.at<double>(2,1) = abs(v[7][s]);
        H.at<double>(0,2) = (-1)*abs(v[2][s]);
        H.at<double>(1,2) = (-1)*abs(v[5][s]);
        H.at<double>(2,2) = abs(v[8][s]);
     
        Mat tmp2=t2.inv();
        Mat tmp1=tmp2*H;
        H=tmp1*t1;
         
        H.at<double>(0,2) = (1)*abs(H.at<double>(0,2));
        H.at<double>(1,2) = (0.5)* abs( H.at<double>(1,2)) ;
        //H.at<double>(3,2) = (-1)* abs( H.at<double>(3,2)) ;
        //cout<<"La matriz H2:"<<endl<<H<<"----";       
        
        /*vector<Point2f> obj_corners(4);
        obj_corners[0] = cvPoint(0,0); obj_corners[1] = cvPoint( img_object.cols, 0 );
        obj_corners[2] = cvPoint( img_object.cols, img_object.rows ); obj_corners[3] = cvPoint( 0, img_object.rows );
        vector<Point2f> scene_corners(4);
        perspectiveTransform( obj_corners, scene_corners, H);
        dibMar(img_matches,img_object,scene_corners);*/
        return H;
}


//  mhomoest: Homography estimation from correcpondence between two    //
//            sets of image points (main program).                     //
//  Programmed by Naoya Ohta and Shimizu Yoshiyuki (1999/2/25)         //
//  Computer Science Department, Gunma University                      //
void maxminCoor(vector<Point2f> ptsMaxMin,vector<Point2f> &scene_corners){
  if(ptsMaxMin[0].y>ptsMaxMin[1].y){  
      scene_corners[0].y=ptsMaxMin[0].y;
      scene_corners[1].y=ptsMaxMin[0].y;
  }else{
      scene_corners[0].y=ptsMaxMin[1].y;
      scene_corners[1].y=ptsMaxMin[1].y;
  }
  if(ptsMaxMin[2].y>ptsMaxMin[3].y){  
      scene_corners[2].y=ptsMaxMin[3].y;
      scene_corners[3].y=ptsMaxMin[3].y;
  }else{
      scene_corners[2].y=ptsMaxMin[2].y;
      scene_corners[3].y=ptsMaxMin[2].y;
  }
  if(ptsMaxMin[0].x>ptsMaxMin[3].x){
      scene_corners[0].x=ptsMaxMin[0].x;
      scene_corners[3].x=ptsMaxMin[0].x;
  }else{
      scene_corners[0].x=ptsMaxMin[3].x;
      scene_corners[3].x=ptsMaxMin[3].x;
  }
  if(ptsMaxMin[1].x>ptsMaxMin[2].x){
      scene_corners[1].x=ptsMaxMin[2].x;
      scene_corners[2].x=ptsMaxMin[2].x;
  }else{
      scene_corners[1].x=ptsMaxMin[1].x;
      scene_corners[2].x=ptsMaxMin[1].x;
  }
      
}
 void calcularPers(Mat img_object,Mat img_scene,Mat img_matches,vector<Point2f> obj,vector<Point2f> scene,Mat H){
 

  //perspectiveTransform( obj_corners, scene_corners, H);
  Mat invH=H.inv();
  Mat ptxp = Mat(3,1, CV_64F, double(0));
   Mat ptx = Mat(3,1, CV_64F, double(0));
   int height=img_scene.rows;
   int width=img_scene.cols;
  
   imshow("scena",img_scene);

   cout<<" ----------------";

   vector<Point2f> obj_corners(4);
  obj_corners[0] = cvPoint(0,0); 
  obj_corners[1] = cvPoint( img_object.cols, 0 );
  obj_corners[2] = cvPoint( img_object.cols, img_object.rows ); 
  obj_corners[3] = cvPoint( 0, img_object.rows );
  vector<Point2f> scene_corners(4);
  vector<Point2f> ptsMaxMin(4);
   
  for(int i=0;i<4;i++){
      ptsMaxMin[i].x=H.at<double>(0,0)*obj_corners[i].x+H.at<double>(0,1)*obj_corners[i].y+H.at<double>(0,2);
      ptsMaxMin[i].y=H.at<double>(1,0)*obj_corners[i].x+H.at<double>(1,1)*obj_corners[i].y+H.at<double>(1,2);
      double val=H.at<double>(2,0)*obj_corners[i].x+H.at<double>(2,1)*obj_corners[i].y+H.at<double>(2,2);
      ptsMaxMin[i].x= ptsMaxMin[i].x/val;
      ptsMaxMin[i].y= ptsMaxMin[i].y/val;
    //  cout<<ptsMaxMin[i].x<<" , "<<ptsMaxMin[i].y<<endl;
  }
  
maxminCoor(ptsMaxMin,scene_corners);
int resFil=int(scene_corners[2].y-scene_corners[0].y);
int resCol=int(scene_corners[2].x-scene_corners[0].x);
 Mat res=Mat::zeros(abs(resFil),abs(resCol),double(0));
   Mat cambio=Mat::zeros(abs(resFil),abs(resCol),double(0));
   for (int i=0; i<img_object.rows; i++){ //y - ver 
        for (int j=0; j<img_object.cols; j++){ //x - hor 
            // set X_a 
            ptxp.at<double>(0,0)=(double)j;
            ptxp.at<double>(1,0)=(double)i;
            ptxp.at<double>(2,0)=1.0;
            // compute X 
            ptx=H*ptxp; 
            //#define CLIP2(minv, maxv, value) (min(maxv, max(minv, value))) 
            int curpi = CLIP2(0, height-1, (int)(ptx.at<double>(1,0)/ptx.at<double>(2,0))); 
            int curpj = CLIP2(0, width-1, (int)(ptx.at<double>(0,0)/ptx.at<double>(2,0))); 

            //cvSet2D(img_object,curpi,curpj,cvGet2D(img_scene,i,j)); 
            //cout<<" "<<H.at<double>(0,2)/H.at<double>(2,2)<<" , "<<H.at<double>(1,2)/H.at<double>(2,2)<<endl;
            if(int((curpi)-scene_corners[0].y<resFil)&&int((curpj)-scene_corners[0].x)<resCol){
                 
               // res.at<uchar>((int)(curpi)-abs(scene_corners[0].y),(int)(curpj)-abs(scene_corners[0].x))=img_object.at<uchar>(i,j);
                res.at<uchar>((int)(curpi)-abs(scene_corners[0].y),(int)(curpj)-abs(scene_corners[0].x))=abs(img_object.at<uchar>(i,j)-img_scene.at<uchar>(curpi,curpj));
            }
        } 
    }


  
  //perspectiveTransform( obj_corners, scene_corners, H);
  
  circle(img_scene,scene_corners[0],5,Scalar(255,255,255),CV_FILLED,8);  
  circle(img_scene,scene_corners[1],4,Scalar(255,255,255),CV_FILLED,8);  
  circle(img_scene,scene_corners[2],4,Scalar(255,255,255),CV_FILLED,8); 
  circle(img_scene,scene_corners[3],4,Scalar(255,255,255),CV_FILLED,8); 
  
/*  circle(res,scene_corners[0],4,Scalar(255,255,255),CV_FILLED,8); 
  circle(res,scene_corners[1],4,Scalar(255,255,255),CV_FILLED,8); 
  circle(res,scene_corners[2],4,Scalar(255,255,255),CV_FILLED,8); 
  circle(res,scene_corners[3],4,Scalar(255,255,255),CV_FILLED,8); */
 // dibMar(img_matches,img_object,scene_corners);
     imshow("imagen de cambio",res);
    imshow( "coordenadas", img_scene );
  waitKey(0);
 }

int main (int argc, char *argv[]){
    cout<<"-----------";
    clock_t start = clock();
    Mat img_object = imread("/home/fredy/imagenes/images7.jpg", CV_LOAD_IMAGE_GRAYSCALE );
    Mat img_scene = imread( "/home/fredy/imagenes/images4.jpg", CV_LOAD_IMAGE_GRAYSCALE );
    vector<Point2f> obj;
    vector<Point2f> scene;
    Mat img_matches=detecDesc(img_object,img_scene,obj, scene);
     Mat H;
    start = clock();
    H=homografia_0(img_object,img_matches,obj,scene);
    cout<<endl<<"El error es: "<<calcularError(H,obj,scene)<<endl;
    calcTime(" 0 ",start);
    calcularPers(img_object,img_scene,img_matches,obj,scene,H);
    /*start = clock();
    H=homografia_RANSAC(img_object,img_matches,obj,scene);
    cout<<endl<<"El error es: "<<calcularError(H,obj,scene)<<endl;
    calcTime(" RANSAC ",start);
     //persTrasn(img_object,img_matches,H);
    
    start = clock();
    H=homografia_LMEDS(img_object,img_matches,obj,scene);
    cout<<endl<<"El error es: "<<calcularError(H,obj,scene)<<endl;
    calcTime(" LMEDS ",start);
     // persTrasn(img_object,img_matches,H);
     start = clock();
    H=dlt2d(img_object,img_matches,obj,scene);
    cout<<endl<<"El error es: "<<calcularError(H,obj,scene)<<endl;
    calcTime(" DLT ",start);
    persTrasn(img_object,img_matches,H);
    //persTrasn(img_object,img_matches,H);
    start = clock();
    H=dlt2DNorm(img_object,img_matches,obj,scene);
    cout<<endl<<"El error es No: "<<calcularError(H,obj,scene)<<endl;
    calcTime(" DLT Norm ",start);
     
    cout<<" "<<kanatani();
    ifstream read("/home/fredy/NetBeansProjects/kanatani/correspSIFTa/results.txt");
    double dato1,dato2,dato3;
    cout<<" //////////////";
    H=Mat(3,3, CV_64F, double(0));
    int i=0;
   while(read>>dato1){
            read>>dato2;
            read>>dato3;
            H.at<double>(i,0)=dato1;
            H.at<double>(i,1)=dato2;
            H.at<double>(i,2)=dato3;
            cout<<" "<<dato1<<" , "<<dato2<<" , "<<dato3<<endl;
            i++;
   }
    //persTrasn(img_object,img_matches,H);
   read.close();*/
    //getchar();
    return 0;
}
 