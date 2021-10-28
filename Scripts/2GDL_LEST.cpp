#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <algorithm>
#include <vector>
#include <opencv2/features2d.hpp>
#include <opencv2/calib3d.hpp>
#include <opencv2/imgproc.hpp>

#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <opencv2/core/eigen.hpp>

#define GLOG_NO_ABBREVIATED_SEVERITIES
#include "ceres/ceres.h"
#include "gflags/gflags.h"
#include "glog/logging.h"



# define M_PI           3.14159265358979323846

#define NIMGS 15

using namespace cv;
using namespace std;


int thresh_max, thresh_val;
double alpha, beta, gamma, tx, ty, tz;
double l1, l2;
double cP2[3], cP3[3];
Mat src, threshold_img, distImage;

Eigen::Matrix3f Triangulation(vector<vector<Point2f>> imPoints, Eigen::MatrixXd c2Tc1, Eigen::MatrixXd c3Tc1);
void AxisDrawing(double alpha, double beta, double gamma, Mat tvec, int index, double q2, vector<Point2f> &imPoints_repro, double L1, double L2);
double norm_dist(Eigen::MatrixXf Point1, Eigen::MatrixXf Point2);

DEFINE_string(minimizer,
	"trust_region",
	"Minimizer type to use, choices are: line_search & trust_region");

//FUNCIONES PARA LA OPTIMIZACIÓN DE MINIMOS CUADRADOS NO LINEALES
struct F1 {
	template <typename T>
	bool operator()(const T* const a, const T* const b, const T* const g, T* residual) const {
		// f1 = L1*(cos(a)*sin(b)*cos(g)+sin(a)*sin(g))+tx-cP2(1);
		residual[0] = l1 * (cos(a[0]) * sin(b[0]) * cos(g[0]) + sin(a[0]) * sin(g[0])) + tx - cP2[0];
		return true;
	}
};

struct F2 {
	template <typename T>
	bool operator()(const T* const a, const T* const b, const T* const g, T* residual) const {
		// f2 = L1*(sin(a)*sin(b)*cos(g)-cos(a)*sin(g))+ty-cP2(2);
		residual[0] = l1 * (sin(a[0]) * sin(b[0]) * cos(g[0]) - cos(a[0]) * sin(g[0])) + ty - cP2[1];
		return true;
	}
};

struct F3 {
	template <typename T>
	bool operator()(const T* const b, const T* const g, T* residual) const {
		// f3 = L1*cos(b)*cos(g)+tz-cP2(3);
		residual[0] = l1 * cos(b[0]) * cos(g[0]) + tz - cP2[2];
		return true;
	}
};

struct F4 {
	template <typename T>
	bool operator()(const T* const a, const T* const b, const T* const g, const T* const q2, T* residual) const {
		// f4 = L2*sin(q2)*cos(a)*cos(b)+(L1+L2*cos(q2))*(cos(a)*sin(b)*cos(g)+sin(a)*sin(g))+tx-cP3(1);
		residual[0] = l2 * sin(q2[0]) * cos(a[0]) * cos(b[0]) + (l1 + l2 * cos(q2[0])) * (cos(a[0]) * sin(b[0]) * cos(g[0]) + sin(a[0]) * sin(g[0])) + tx - cP3[0];
		return true;
	}
};

struct F5 {
	template <typename T>
	bool operator()(const T* const a, const T* const b, const T* const g, const T* const q2, T* residual) const {
		// f5 = L2*sin(q2)*sin(a)*cos(b)+(L1+L2*cos(q2))*(sin(a)*sin(b)*cos(g)-cos(a)*sin(g))+ty-cP3(2);
		residual[0] = l2 * sin(q2[0]) * sin(a[0]) * cos(b[0]) + (l1 + l2 * cos(q2[0])) * (sin(a[0]) * sin(b[0]) * cos(g[0]) - cos(a[0]) * sin(g[0])) + ty - cP3[1];
		return true;
	}
};

struct F6 {
	template <typename T>
	bool operator()(const T* const b, const T* const g, const T* const q2, T* residual) const {
		// f6 = L2*sin(q2)*(-sin(b))+(L1+L2*cos(q2))*cos(b)*cos(g)+tz-cP3(3)]);
		residual[0] = l2 * sin(q2[0]) * (-sin(b[0])) + (l1 + l2 * cos(q2[0])) * cos(b[0]) * cos(g[0]) + tz - cP3[2];
		return true;
	}
};

int main(int argc, char** argv)
{
	GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
	google::InitGoogleLogging(argv[0]);

	ofstream myfile;
	myfile.open("Total_time_2gdl.txt", ios::out | ios::trunc);
	myfile.close();
	myfile.open("nIterations_2gdl.txt", ios::out | ios::trunc);
	myfile.close();
	myfile.open("Full_report2gdl.txt", ios::out | ios::trunc);
	myfile.close();
	myfile.open("P3D_2gdl.txt", ios::out | ios::trunc);
	myfile.close();
	myfile.open("2gdl_repro.txt", ios::out | ios::trunc);
	myfile.close();

	// Transformaciones
	Eigen::MatrixXd c2Tc1(3, 4);
	c2Tc1 << 0.1696, 0.5295, -0.8312, 2.2852,
		-0.6679, 0.6819, 0.2981, -0.8720,
		 0.7247, 0.5046, 0.4693, 1.5080;

	Eigen::MatrixXd c3Tc1(3, 4);
	c3Tc1 << -0.2186, -0.8543, 0.4716, -1.2109,
		 0.5033, 0.3154, 0.8045, -2.1316,
		-0.8360, 0.4132, 0.3610, 1.4416;

	for (int k = 1; k <= NIMGS+1; k++) { //Bucle principal

		//Cargar imagenes -> si orden cambia, cambiar aqui
		char path1[99]; char path2[99]; char path3[99];
		sprintf_s(path1, 99, "2gdl/2-%d.bmp", k);
		sprintf_s(path2, 99, "2gdl/1-%d.bmp", k);
		sprintf_s(path3, 99, "2gdl/3-%d.bmp", k);

		printf("\nIMAGENES - %d\n", k);
		// Load the image
		string path;
		string paths[] = { path1,path2,path3 };
		char aux[99];
		int i;
		int crop_top[] = { 300,280,200 };
		int crop_bot[] = { 0,40,120 };
		int crop_left[] = { 100,100,100 };
		int crop_right[] = { 230,170,170 };
		int im_validas = 0;
		vector<vector<Point2f>> imPoints(3);

		Eigen::Matrix4d c1Tw_promed;

		//Scalar color(rand() & 255, rand() & 255, rand() & 255);
		vector<Scalar> colores{ Scalar(255,0,0), Scalar(0,255,0), Scalar(0,0,255) };

		for (i = 0; i < 3; i++) {
			path = paths[i];
			Mat src = imread(path);
			if (src.empty())
			{
				return -1;
			}

			const Rect roi(crop_left[i], crop_top[i], src.cols - crop_right[i] - crop_left[i], src.rows - crop_bot[i] - crop_top[i]);
			src = src(roi).clone();
			Mat original = src.clone();

			Mat mask;
			inRange(src, Scalar(255, 255, 255), Scalar(255, 255, 255), mask);
			threshold(mask, mask, 180, 255, THRESH_BINARY_INV);
			src.setTo(Scalar(0, 0, 0), mask);

			Mat kernel = (Mat_<float>(3, 3) <<
				1, 1, 1,
				1, -8, 1,
				1, 1, 1);

			Mat imgLaplacian;
			filter2D(src, imgLaplacian, CV_32F, kernel);
			Mat sharp;
			src.convertTo(sharp, CV_32F);
			Mat imgResult = sharp - imgLaplacian;
			// convert back to 8bits gray scale
			imgResult.convertTo(imgResult, CV_8UC3);
			cvtColor(imgResult, imgResult, COLOR_BGR2GRAY);

			threshold(imgResult, imgResult, 240, 255, THRESH_BINARY | THRESH_OTSU);

			Mat transfMorfolog = imgResult.clone();
			Mat kernel1 = Mat::ones(5, 5, CV_8U);
			Mat kernel2 = Mat::ones(3, 3, CV_8U);

			//erode(transfMorfolog, transfMorfolog, kernel1);
			dilate(transfMorfolog, transfMorfolog, kernel1);

			Mat transfMorfolog_8u;
			transfMorfolog.convertTo(transfMorfolog_8u, CV_8U);
			// Find total markers
			vector<vector<Point>> contours;
			vector<vector<Point>> contours_circ;
			findContours(transfMorfolog_8u, contours_circ, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);

			//Filtrado por circularidad
			if (contours_circ.size() > 3) {
				vector<double> circularity(contours_circ.size());
				vector<double> circularity_sorted(contours_circ.size());

				for (int j = 0; j < contours_circ.size(); j++) {
					double perimeter = arcLength(contours_circ[j], true);
					double area = contourArea(contours_circ[j]);

					if (perimeter > 0) {
						circularity[j] = 4 * M_PI * (area / (perimeter * perimeter));
					}
				}

				circularity_sorted = circularity;
				sort(circularity_sorted.begin(), circularity_sorted.end(), greater<double>());
				for (int j = 0; j < 3; j++) {
					int s = find(circularity.begin(), circularity.end(), circularity_sorted[j]) - circularity.begin();
					contours.push_back(contours_circ[s]);
				}
			}
			else
			{
				contours = contours_circ;
			}

			// Create the marker image
			Mat markers = Mat::zeros(transfMorfolog.size(), CV_32S);
			// Draw the markers
			for (size_t j = 0; j < contours.size(); j++)
			{
				drawContours(markers, contours, static_cast<int>(j), Scalar(static_cast<int>(j) + 1), -1);
			}
			// Draw the background marker
			circle(markers, Point(5, 5), 3, Scalar(255), -1);
			Mat markers8u;
			markers.convertTo(markers8u, CV_8U, 10);

			int ncomp = contours.size();
			Mat drawing = original.clone();
			if (ncomp > 0) {
				vector<Moments> mu(ncomp);
				vector<Point2f> mc(ncomp); //Mass centers
				vector<Point2f> mc_ordered(ncomp); //Ordered mass centers
				vector<float> area(ncomp);
				vector<float> pos_X(ncomp);
				vector<float> pos_Y(ncomp);
				vector<vector<Point>> finalContours;

				for (int seed = 1; seed <= contours.size(); ++seed)
				{
					cv::Mat1b mask = (markers == seed);

					findContours(mask, finalContours, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);
					mu[seed - 1] = moments(finalContours[0], false);
					mc[seed - 1] = Point2f(mu[seed - 1].m10 / mu[seed - 1].m00, mu[seed - 1].m01 / mu[seed - 1].m00); //Centro de masas
					pos_X[seed - 1] = mc[seed - 1].x;
					pos_Y[seed - 1] = mc[seed - 1].y;
					area[seed - 1] = mu[seed - 1].m00; //Area

					Point2f centro; float radio;
					minEnclosingCircle(finalContours[0], centro, radio);

					pos_X[seed - 1] = centro.x;
					pos_Y[seed - 1] = centro.y;
					circle(drawing, centro, radio, Scalar(255, 0, 255));
					drawContours(drawing, finalContours, -1, Scalar(120, 120, 120), 1, 8, noArray(), 0, Point());
				}

				if (ncomp == 3) {
					int ord_vec[] = { max_element(area.begin(),area.end()) - area.begin(), -1, min_element(area.begin(),area.end()) - area.begin() }; //Ordenamos segun el area
					vector<int> indexes = { 0,1,2 };
					remove_copy(indexes.begin(), indexes.end(), indexes.begin(), ord_vec[0]);
					remove_copy(indexes.begin(), indexes.end(), indexes.begin(), ord_vec[2]);
					ord_vec[1] = indexes[0];

					for (int j = 0; j < contours.size(); j++) {
						Scalar color = colores[j];
						circle(drawing, mc[ord_vec[j]], 4, color, -1, 8, 0);
						imPoints[i].push_back(Point2f(pos_X[ord_vec[j]] + crop_left[i], pos_Y[ord_vec[j]] + crop_top[i]));																													   //sumar los recortes en px
					}
					im_validas++;
					line(drawing, Point2f(pos_X[ord_vec[0]],pos_Y[ord_vec[0]]), Point2f(pos_X[ord_vec[1]], pos_Y[ord_vec[1]]), Scalar(0, 255, 255), 1);
					line(drawing, Point2f(pos_X[ord_vec[1]], pos_Y[ord_vec[1]]), Point2f(pos_X[ord_vec[2]], pos_Y[ord_vec[2]]), Scalar(0, 255, 255), 1);
				}
			}
			sprintf_s(aux, 99, "Final - %d", i + 1);
			imshow(aux, drawing);
			printf("\nncomp = %d\tImage %d\n", ncomp, i);
		}


		if (im_validas == 3) { //Si las 3 imagenes son validas --> Triangulamos
			cout << "\n\n VALID IMAGE SET" << endl;

			Eigen::Matrix3f points3D;

			for (int i = 0; i < 3; i++) {
				sprintf_s(aux, 99, "centers_resultsC%d.txt", i + 1);
				myfile.open(aux, ios::out | ios::app);
				for (int ix = 0; ix < 3; ix++) {
					myfile << imPoints[i][ix].x << "," << imPoints[i][ix].y << endl;
				}

				myfile.close();
			}

			//Hacemos la estimación de los puntos 3D mediante triangulación, vistos desde la camara 1
			points3D = Triangulation(imPoints,c2Tc1,c3Tc1);

			//Resolvemos mediante mínimos cuadradados cP = cTobj * objP
			//Damos a tx, ty y tz los puntos del primer marcador para simplificar el problema:
			tx = points3D(0,0);
			ty = points3D(1,0);
			tz = points3D(2,0);

			//Cargamos los valores de los puntos de los demás marcadores:
			cP2[0] = points3D(0,1);
			cP2[1] = points3D(1,1);
			cP2[2] = points3D(2,1);

			cP3[0] = points3D(0,2);
			cP3[1] = points3D(1,2);
			cP3[2] = points3D(2,2);

			//Cargamos las estimaciones de las longitudes de las articulaciones
			l1 = norm_dist(points3D.block(0, 0, 3, 1), points3D.block(0, 1, 3, 1));
			l2 = norm_dist(points3D.block(0, 1, 3, 1), points3D.block(0, 2, 3, 1));

			//Damos las condiciones inciales para las variables
			double alpha = 0.0;
			double beta = 0.0;
			double gamma = 0.0;
			double q2 = 0.0;

			ceres::Problem problem;
			problem.AddResidualBlock(new ceres::AutoDiffCostFunction<F1, 1, 1, 1, 1>(new F1), NULL, &alpha, &beta, &gamma);
			problem.AddResidualBlock(new ceres::AutoDiffCostFunction<F2, 1, 1, 1, 1>(new F2), NULL, &alpha, &beta, &gamma);
			problem.AddResidualBlock(new ceres::AutoDiffCostFunction<F3, 1, 1, 1>(new F3), NULL, &beta, &gamma);
			problem.AddResidualBlock(new ceres::AutoDiffCostFunction<F4, 1, 1, 1, 1, 1>(new F4), NULL, &alpha, &beta, &gamma, &q2);
			problem.AddResidualBlock(new ceres::AutoDiffCostFunction<F5, 1, 1, 1, 1, 1>(new F5), NULL, &alpha, &beta, &gamma, &q2);
			problem.AddResidualBlock(new ceres::AutoDiffCostFunction<F6, 1, 1, 1, 1>(new F6), NULL, &beta, &gamma, &q2);
			
			problem.SetParameterLowerBound(&alpha, 0, -2 * M_PI);
			problem.SetParameterLowerBound(&beta, 0, -2 * M_PI);
			problem.SetParameterLowerBound(&gamma, 0, -2 * M_PI);
			problem.SetParameterLowerBound(&q2, 0, -M_PI);

			problem.SetParameterUpperBound(&alpha, 0, 2 * M_PI);
			problem.SetParameterUpperBound(&beta, 0, 2 * M_PI);
			problem.SetParameterUpperBound(&gamma, 0, 2 * M_PI);
			problem.SetParameterUpperBound(&q2, 0, M_PI);

			ceres::Solver::Options options;
			LOG_IF(
				FATAL,
				!ceres::StringToMinimizerType(FLAGS_minimizer, &options.minimizer_type))
				<< "Invalid minimizer: " << FLAGS_minimizer
				<< ", valid options are: trust_region and line_search.";

			options.max_num_iterations = 50;
			options.linear_solver_type = ceres::DENSE_QR;
			options.minimizer_progress_to_stdout = true;

			// Ejecutamos el Solver
			ceres::Solver::Summary summary;
			ceres::Solve(options, &problem, &summary);

			cout << summary.FullReport() << "\n";

			// clang-format off
			cout << "Final alpha = " << alpha * 180 / M_PI
				<< ", beta = " << beta * 180 / M_PI
				<< ", gamma = " << gamma * 180 / M_PI
				<< ", q2 = " << q2 * 180 / M_PI
				<< ", l1 = " << l1
				<< ", l2 = " << l2
				<< "\n";
			
			vector<Point2f> imPoints_reproyected;
			Mat tvec = (Mat_<double>(3, 1) << tx, ty, tz);
			AxisDrawing(alpha,beta,gamma,tvec,k,q2,imPoints_reproyected,l1,l2);

			sprintf_s(aux, 99, "2gdl_repro.txt");
			myfile.open(aux, ios::out | ios::app);
			for (int i = 0; i < 3; i++) {
				myfile << imPoints_reproyected[i].x << "," << imPoints_reproyected[i].y << endl;
			}
			myfile.close();

			sprintf_s(aux, 99, "Total_time_2gdl.txt");
			myfile.open(aux, ios::out | ios::app);
			myfile << summary.total_time_in_seconds << endl;
			myfile.close();

			sprintf_s(aux, 99, "nIterations_2gdl.txt");
			myfile.open(aux, ios::out | ios::app);
			myfile << summary.num_successful_steps + summary.num_unsuccessful_steps << endl;
			myfile.close();

			sprintf_s(aux, 99, "Full_report2gdl.txt");
			myfile.open(aux, ios::out | ios::app);
			myfile << alpha * 180 / M_PI << " " << beta * 180 / M_PI << " " << gamma * 180 / M_PI << " " << q2 * 180 / M_PI << endl;
			myfile.close();

			cout << "Press any key to continue" << endl;
			waitKey(0);
		}
	}
	return 0;
}


Eigen::Matrix3f Triangulation(vector<vector<Point2f>> imPoints, Eigen::MatrixXd c2Tc1, Eigen::MatrixXd c3Tc1)
{
	ofstream myfile; char aux[99];

	Mat cameraMatrix1 = (Mat_<float>(3, 3) <<
		943.3056, 0, 681.2714,
		0, 947.5338, 489.8241,
		0, 0, 1.0000);
	Mat distCoeffs1 = (Mat_<float>(4, 1) << -0.2713, 0.0183, 0, 0);

	Mat cameraMatrix2 = (Mat_<float>(3, 3) <<
		965.5391, 0, 613.7441,
		0, 965.8018, 519.5804,
		0, 0, 1.0000);
	Mat distCoeffs2 = (Mat_<float>(4, 1) << -0.3045, 0.1392, 0, 0);

	Mat cameraMatrix3 = (Mat_<float>(3, 3) <<
		959.0305, 0, 678.4766,
		0, 956.9934, 518.3325,
		0, 0, 1.0000);
	Mat distCoeffs3 = (Mat_<float>(4, 1) << -0.3149, 0.1334, 0, 0);
	Mat projMatrC1 = (Mat_<float>(3, 4) <<
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0); //Primera matriz de proyeccion para la camara fija 1, que debe ser nula

	cameraMatrix2.convertTo(cameraMatrix2, CV_64F); //Hacemos un cambio del tipo de dato para que puedan
	cameraMatrix3.convertTo(cameraMatrix3, CV_64F); //ser multiplicadas

	Mat points4D_c1c2, points4D_c1c3;
	Mat projMatrC2, projMatrC3;
	eigen2cv(c2Tc1, projMatrC2);
	eigen2cv(c3Tc1, projMatrC3); //Matrices de proyeccion para las camaras 2 y 3

	triangulatePoints(cameraMatrix1 * projMatrC1, cameraMatrix2 * projMatrC2, imPoints[0], imPoints[1], points4D_c1c2);
	triangulatePoints(cameraMatrix1 * projMatrC1, cameraMatrix3 * projMatrC3, imPoints[0], imPoints[2], points4D_c1c3);

	Eigen::MatrixXf points3D_c1c2, points3D_c1c3, points3D;
	cv2eigen(points4D_c1c2, points3D_c1c2); //Pasamos a Eigen, para que sea mas sencillo operar
	cv2eigen(points4D_c1c3, points3D_c1c3);

	for (int k = 0; k < 3; k++) {
		points3D_c1c2.block(0, k, 4, 1) = points3D_c1c2.block(0, k, 4, 1) / points3D_c1c2(3, k);
		points3D_c1c3.block(0, k, 4, 1) = points3D_c1c3.block(0, k, 4, 1) / points3D_c1c3(3, k);
	}
	points3D = (points3D_c1c2.block(0,0,3,3) + points3D_c1c3.block(0,0,3,3)) / 2.0;

	sprintf_s(aux, 99, "P3D_2gdl.txt"); //Guardamos en un archivo para estudiar en MATLAB
	myfile.open(aux, ios::out | ios::app);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			myfile << points3D_c1c2(j, i) << "," << points3D_c1c3(j, i) << endl;
		}
	}
	myfile.close();
	cout << "Points3D: " << points3D << endl;

	return points3D;
}

void AxisDrawing(double alpha, double beta, double gamma, Mat tvec, int index, double q2, vector<Point2f>&imPoints_repro, double L1, double L2)
{
	Eigen::Matrix3d Rotx, Roty, Rotz, Rot;
	Mat rvec;
	Mat cameraMatrix = (Mat_<float>(3, 3) <<
		943.3056, 0, 681.2714,
		0, 947.5338, 489.8241,
		0, 0, 1.0000);
	Mat distCoeffs = (Mat_<float>(4, 1) << -0.2713, 0.0183, 0, 0);

	vector<Point3f> object_AxisPoints = { Point3f(0.5,0,0), Point3f(0,0.5,0), Point3f(0,0,0.5), Point3f(0,0,0) };
	vector<Point3f> objectPoints_repro = { Point3f(0,0,0), Point3f(0,0,L1), Point3f(L2 * sin(q2),0,L1 + L2 * cos(q2)) };
	vector<Point2f> image_AxisPoints;

	Rotz << cos(alpha), -sin(alpha), 0,
			sin(alpha), cos(alpha), 0,
			0, 0, 1;
	Roty << cos(beta), 0, sin(beta),
			0, 1, 0,
			-sin(beta), 0, cos(beta);
	Rotx << 1, 0, 0,
		0, cos(gamma), -sin(gamma),
		0, sin(gamma), cos(gamma);

	Rot = Rotz * Roty * Rotx;
	eigen2cv(Rot, rvec);
	Rodrigues(rvec, rvec);

	projectPoints(objectPoints_repro, rvec, tvec, cameraMatrix, distCoeffs, imPoints_repro);
	projectPoints(object_AxisPoints, rvec, tvec, cameraMatrix, distCoeffs, image_AxisPoints);

	char path[99];
	sprintf_s(path, 99, "2gdl/2-%d.bmp", index);
	Mat src = imread(path);

	//Dibujamos la lineas
	line(src, image_AxisPoints[3], image_AxisPoints[0], Scalar(0, 0, 255), 2); //Eje X rojo
	line(src, image_AxisPoints[3], image_AxisPoints[1], Scalar(0, 255, 0), 4); //Eje Y verde
	line(src, image_AxisPoints[3], image_AxisPoints[2], Scalar(255, 0, 0), 6); //Eje Z azul
	char line[99];
	sprintf_s(line, 99, "q2=%f", q2*180/M_PI);
	putText(src,line,Point(100,100),FONT_HERSHEY_COMPLEX_SMALL,1,Scalar(255, 255, 255), 1,false);

	imshow("Axis", src);
}

double norm_dist(Eigen::MatrixXf Point1, Eigen::MatrixXf Point2) {
	double dist;
	dist = sqrt(pow(Point1(0, 0) - Point2(0, 0), 2) + pow(Point1(1, 0) - Point2(1, 0), 2) + pow(Point1(2, 0) - Point2(2, 0), 2));

	return dist;
}