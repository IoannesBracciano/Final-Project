#pragma once

/***
 * FBLF.HPP
 * Fast Bilateral Filter
 * A fast approximation of the bilateral filter using 3D convolution and a truncated kernel.
 * This implementation is strongly based on the paper of Paris and Durand, 2009.
 *
 * Author: Ioannes Bracciano <john.bracciano@gmail.com>
 * All rights reserved
 */

#include <iostream>

 // OpenCV libraries
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>

#define DEBUG 0
#define SUPPORTED_CHANNELS 1
#define KERNEL_SIZE 5
#define KERNEL_SIGMA 1

#define DIM_Y 1;
#define DIM_X 2;
#define DIM_Z 3;

using namespace cv;

template <typename _Tp, int n>
Vec<_Tp, n> round(const Vec<_Tp, n> vec) {
	Vec<_Tp, n> roundvec;

	for (int i = 0; i < n; ++i) {
		roundvec[i] = round(vec[i]);
	}

	return roundvec;
}
//
//void pad(int padding, InputArray &in, uchar dim) {
//	Vec<int, 3> padding_3rd = { 0, 0, padding };
//	Range ranges_3rd[3] = {
//		Range::all(), Range::all(), Range(padding, padding + range_)
//	};
//	Mat w_i_pad_3rd = Mat(3, (int*)Mat(padding_3rd + sizes + padding_3rd).data, CV_16U, Scalar(0));
//}
//
//void convolve(const Vec<double, KERNEL_SIZE> &kernel, InputArray in, OutputArray out, uchar dim) {
//
//}

double interpolate(Mat V, double y, double x, double z) {
	static int height_ = V.size[0];
	static int width_ = V.size[1];
	static int range_ = V.size[2];
	static int step_y = width_ * range_;
	static int step_x = range_;

	if (height_ != V.size[0] || width_ != V.size[1] || range_ != V.size[2]) {
		height_ = V.size[0];
		width_ = V.size[1];
		range_ = V.size[2];
		step_y = width_ * range_;
		step_x = range_;
	}

	/*int height_ = V.size[0];
	int width_ = V.size[1];
	int range_ = V.size[2];
	int step_y = width_ * range_;
	int step_x = range_;*/

	/*std::cout << "{height_, width_, range_}: {" << height_ << ", " << width_ << ", " << range_ << "}" << std::endl;
	std::cout << "{qy, qx, qz}: {" << y << ", " << x << ", " << z << "}" << std::endl;
	getchar();*/

	int y0 = floor(y);    int y1 = y0 + 1;    double yd = (y - y0);
	int x0 = floor(x);    int x1 = x0 + 1;    double xd = (x - x0);
	int z0 = floor(z);    int z1 = z0 + 1;    double zd = (z - z0);

	int idx000[3] = { y0, x0, z0 };    int idx010[3] = { y0, x1, z0 };
	int idx100[3] = { y1, x0, z0 };    int idx110[3] = { y1, x1, z0 };
	int idx001[3] = { y0, x0, z1 };    int idx011[3] = { y0, x1, z1 };
	int idx101[3] = { y1, x0, z1 };    int idx111[3] = { y1, x1, z1 };

	int y0_linear = y0 * step_y;    int y1_linear = y1 * step_y;
	int x0_linear = x0 * step_x;    int x1_linear = x1 * step_x;

	int idx000_linear = y0_linear + x0_linear + z0;
	int idx010_linear = y0_linear + x1_linear + z0;
	int idx100_linear = y1_linear + x0_linear + z0;
	int idx110_linear = y1_linear + x1_linear + z0;
	int idx001_linear = y0_linear + x0_linear + z1;
	int idx011_linear = y0_linear + x1_linear + z1;
	int idx101_linear = y1_linear + x0_linear + z1;
	int idx111_linear = y1_linear + x1_linear + z1;

	double *V_p = (double*)V.data;

	double V000 = V.isContinuous() ? *(V_p + idx000_linear) : V.at<double>(idx000);
	double V010 = V.isContinuous() ? *(V_p + idx010_linear) : V.at<double>(idx010);
	double V100 = V.isContinuous() ? *(V_p + idx100_linear) : V.at<double>(idx100);
	double V110 = V.isContinuous() ? *(V_p + idx110_linear) : V.at<double>(idx110);
	double V001 = V.isContinuous() ? *(V_p + idx001_linear) : V.at<double>(idx001);
	double V011 = V.isContinuous() ? *(V_p + idx011_linear) : V.at<double>(idx011);
	double V101 = V.isContinuous() ? *(V_p + idx101_linear) : V.at<double>(idx101);
	double V111 = V.isContinuous() ? *(V_p + idx111_linear) : V.at<double>(idx111);

	/*if(y==0.5)
	std::cout << "{ "	<< V000 << "," << V010 << ","
	<< V100 << "," << V110 << ","
	<< V001 << "," << V011 << ","
	<< V101 << "," << V111 << " }" << std::endl;*/

	double c00 = V000*(1 - xd) + V010*xd;
	double c10 = V100*(1 - xd) + V110*xd;
	double c01 = V001*(1 - xd) + V011*xd;
	double c11 = V101*(1 - xd) + V111*xd;

	double c0 = c00*(1 - yd) + c10*yd;
	double c1 = c01*(1 - yd) + c11*yd;

	double c = c0*(1 - zd) + c1*zd;

	return c;
}

static void extractZSliceOBL(const cv::Mat& image3d, const int z, cv::Mat& slice)
{
	// create the roi 
	cv::Range ranges[3];
	ranges[0] = cv::Range::all();
	ranges[1] = cv::Range::all();
	ranges[2] = cv::Range(z, z + 1);

	// get the roi from the image; 
	// calling clone() makes sure the data is continuous 
	slice = image3d(ranges).clone();

	// create a temporarily 2d image and copy its size 
	// to make our image a real 2d image 
	cv::Mat temp2d;
	temp2d.create(2, &(image3d.size[0]), image3d.type());
	slice.copySize(temp2d);
}


double fastBilateralFilter(InputArray _src, OutputArray _dst, float sigmaSpatial, float sigmaRange, int intensities = 256) {
	Mat src = _src.getMat();

	uchar * image_p = src.data;
		// Pointer to the raw data of the image
	CV_Assert(image_p != NULL);
		// Break if image could not be loaded

	const int height = src.rows;
	const int width = src.cols;
	const int range = intensities;

	const int height_ = round(height / sigmaSpatial) + 1;
	const int width_ = round(width / sigmaSpatial) + 1;
	const int range_ = round(range / sigmaRange) + 1;

	const uint channels = src.channels();
	// If image has less or more than 3 channels, stop and show error
	CV_Assert(channels == SUPPORTED_CHANNELS);
	// Break if image has more or less channels than SUPPORTED_CHANNELS

	double t = (double)getTickCount();
	// Time is ticking...!

	/*
	* Beginning of the basic algorithm
	* Algorithm steps written exactly as introduced in page 12 of Paris and Durand paper
	*
	* Step 1: Initialize all w_i_ and w_ values to 0
	*/
	if (DEBUG) std::cout << "DEBUG: Initialization..." << std::endl;

	Vec<int, 3> sizes = { height_, width_, range_ };
	Mat w_i_ = Mat(3, (int*)Mat(sizes).data, CV_64F, Scalar(0));
	Mat w_ = Mat(3, (int*)Mat(sizes).data, CV_64F, Scalar(0));

	// Retrieving pointers to the data of the matrices for efficient element access
	// The matrices need to be continues in the memory for that
	double *w_i_p = (double*)w_i_.data;
	double *w_p = (double*)w_.data;
	/*
	* Step 2:  Compute the minimum intensity value
	*/
	if (DEBUG) std::cout << "DEBUG: Computing image minimum..." << std::endl;

	double _min, _max;
	// OpenCV won't let us have the minimum values
	// without having the maximum as well...
	uchar min;

	minMaxLoc(src, &_min, &_max);
	// Getting the min and the max values of each channel
	min = (uchar)_min;
	// Converting the values to unsigned char

	/*
	* Step 3: For each pixel (Y,X) belonging to S with an intensity I(Y,X) belonging to R...
	*/
	if (DEBUG) {
		if (w_i_.isContinuous() && w_.isContinuous() && src.isContinuous()) {
			std::cout << "DEBUG: Downsampling with continuous matrices" << std::endl;
		}
		else {
			std::cout << "DEBUG: Downsampling" << std::endl;
		}
	}

	// Precalculating steps for finding linear indexes
	const unsigned long step_y = width_*range_;
	const unsigned long step_x = range_;
	const unsigned long step_z = 1;

	for (uint Y = 0; Y < height; ++Y) {
		for (uint X = 0; X < width; ++X) {
			// a) Compute the homogeneous vector (wi,w)
			uchar wi;
			if (src.isContinuous()) wi = image_p[Y*width + X];
			else wi = src.at<uchar>(Y, X);
			// b) Compute the downsampled coordinates
			uint y = round(Y / sigmaSpatial);
			uint x = round(X / sigmaSpatial);
			uchar zeta = round((wi) / sigmaRange);
			// c) Update the downsampled SxR space
			int idx[3] = { y, x, zeta };
			int linear_idx = y*step_y + x*step_x + zeta*step_z;
			//std::cout << "Sizes: " << w_i_.size[0] << ", " << w_i_.size[1] << ", " << w_i_.size[2] << ", " << w_i_.size[3] << ", " << w_i_.size[4] << std::endl;
			//std::cout << "Accessing: " << idx[0] << ", " << idx[1] << ", " << idx[2] << std::endl;
			//std::cout << "Lnear index: " << linear_idx << std::endl;
			if (w_i_.isContinuous()) w_i_p[linear_idx] += static_cast<double>(wi);
			else w_i_.at<double>(idx) += static_cast<double>(wi);
			if (w_.isContinuous())	w_p[linear_idx] += static_cast<double>(1);
			else w_.at<double>(idx) += static_cast<double>(1);

			/*std::cout << w_i_.at<UINT16>(idx) << ','
			<< w_.at<UINT16>(idx)
			<< std::endl;
			getchar();*/
		}
	}

	/*for (uchar i = 0; i < range_; ++i) {
		Mat slice;
		extractZSliceOBL(w_i_.mul(w_i_), i, slice);
		imshow("w_i_", slice);
		waitKey(0);
	}*/

	/*
	* Step 3:  Convolve (w_i_, w_) with a 3D Gaussian
	*/
	Vec<double, KERNEL_SIZE> kernel = getGaussianKernel(KERNEL_SIZE, KERNEL_SIGMA);
	const int padding = KERNEL_SIZE / 2;

	Mat w_i_b = Mat(3, (int*)Mat(sizes).data, CV_64F, Scalar(0));
	Mat w_b = Mat(3, (int*)Mat(sizes).data, CV_64F, Scalar(0));
	// Retrieving pointers to data of matrices for quick element access,
	// provided that they are continuous
	double * w_i_b_p = (double*)w_i_b.data;
	double * w_b_p = (double*)w_b.data;

	// Convolving along the 3th dimension
	if (DEBUG) std::cout << "DEBUG: Convolving along the 3th dimension ";
	// Padding image for the convolution
	Vec<int, 3> padding_3rd = { 0, 0, padding };
	Range ranges_3rd[3] = {
		Range::all(), Range::all(), Range(padding, padding + range_)
	};
	Mat w_i_pad_3rd = Mat(3, (int*)Mat(padding_3rd + sizes + padding_3rd).data, CV_64F, Scalar(0));
	Mat w_pad_3rd = Mat(3, (int*)Mat(padding_3rd + sizes + padding_3rd).data, CV_64F, Scalar(0));
	double *w_i_pad_3rd_p = (double*)w_i_pad_3rd.data;
	double *w_pad_3rd_p = (double*)w_pad_3rd.data;

	w_i_.copyTo(w_i_pad_3rd(ranges_3rd));
	w_.copyTo(w_pad_3rd(ranges_3rd));

	// Deleting w_i_ and w_ to free memory
	w_i_.release();
	w_.release();

	// Precalculating linear index steps for each new padded dimension
	const unsigned int step_y_pad_3rd = width_ * (range_ + padding * 2);
	const unsigned int step_x_pad_3rd = (range_ + padding * 2);
	const unsigned int step_z_pad_3rd = 1; // Unchanged

	if (DEBUG) {
		if (w_i_pad_3rd.isContinuous() && w_pad_3rd.isContinuous() &&
			w_i_b.isContinuous() && w_b.isContinuous()) {
			std::cout << "with continuous matrices" << std::endl;
		}
		else {
			std::cout << std::endl;
		}
	}

	for (int y = 0; y < height_; ++y) {
		for (int x = 0; x < width_; ++x) {
			for (int z = 0; z < range_; ++z) {
				Vec<double, KERNEL_SIZE> w_i_elems, w_elems, w_i_elems2;

				if (w_i_pad_3rd.isContinuous() && w_i_b.isContinuous() && w_b.isContinuous()) {
					int linear_idx_pad = y*step_y_pad_3rd + x*step_x_pad_3rd + z*step_z_pad_3rd;
					int linear_idx = y*step_y + x*step_x + z*step_z;

					w_i_elems = {
						w_i_pad_3rd_p[linear_idx_pad + 0 * step_z_pad_3rd],
						w_i_pad_3rd_p[linear_idx_pad + 1 * step_z_pad_3rd],
						w_i_pad_3rd_p[linear_idx_pad + 2 * step_z_pad_3rd],
						w_i_pad_3rd_p[linear_idx_pad + 3 * step_z_pad_3rd],
						w_i_pad_3rd_p[linear_idx_pad + 4 * step_z_pad_3rd]
					};

					w_i_b_p[linear_idx] = kernel.dot(w_i_elems);

					w_elems = {
						w_pad_3rd_p[linear_idx_pad + 0 * step_z_pad_3rd],
						w_pad_3rd_p[linear_idx_pad + 1 * step_z_pad_3rd],
						w_pad_3rd_p[linear_idx_pad + 2 * step_z_pad_3rd],
						w_pad_3rd_p[linear_idx_pad + 3 * step_z_pad_3rd],
						w_pad_3rd_p[linear_idx_pad + 4 * step_z_pad_3rd]
					};

					w_b_p[linear_idx] = kernel.dot(w_elems);

				}
				else {
					int idx1[3] = { y, x, z + 0 };
					int idx2[3] = { y, x, z + 1 };
					int idx3[3] = { y, x, z + 2 };
					int idx4[3] = { y, x, z + 3 };
					int idx5[3] = { y, x, z + 4 };

					w_i_elems = {
						w_i_pad_3rd.at<double>(idx1),
						w_i_pad_3rd.at<double>(idx2),
						w_i_pad_3rd.at<double>(idx3),
						w_i_pad_3rd.at<double>(idx4),
						w_i_pad_3rd.at<double>(idx5)
					};

					w_i_b.at<double>(idx1) = kernel.dot(w_i_elems);

					w_elems = {
						w_pad_3rd.at<double>(idx1),
						w_pad_3rd.at<double>(idx2),
						w_pad_3rd.at<double>(idx3),
						w_pad_3rd.at<double>(idx4),
						w_pad_3rd.at<double>(idx5)
					};

					w_b.at<double>(idx1) = kernel.dot(w_elems);
				} // if - else
			} // 3rd dimension
		} // 2nd dimension
	} // 1st dimension

	  /*for (uchar i = 0; i < range_; ++i) {
		  Mat slice;
		  extractZSliceOBL(w_i_b.mul(.0001), i, slice);
		  imshow("w_i_b_1st", slice);
		  waitKey(0);
	  }*/

	  // Deleting w_i_pad_3rd and w_pad_3rd to free memory
	w_i_pad_3rd.release();
	w_pad_3rd.release();

	// Convolving along the 2th dimension
	if (DEBUG) std::cout << "DEBUG: Convolving along the 2nd dimension ";
	// Padding image for the convolution
	Vec<int, 3> padding_2nd = { 0, padding, 0 };
	Range ranges_2nd[3] = {
		Range::all(), Range(padding, padding + width_), Range::all()
	};
	Mat w_i_pad_2nd = Mat(3, (int*)Mat(padding_2nd + sizes + padding_2nd).data, CV_64F, Scalar(0));
	Mat w_pad_2nd = Mat(3, (int*)Mat(padding_2nd + sizes + padding_2nd).data, CV_64F, Scalar(0));
	double *w_i_pad_2nd_p = (double*)w_i_pad_2nd.data;
	double *w_pad_2nd_p = (double*)w_pad_2nd.data;

	w_i_b.copyTo(w_i_pad_2nd(ranges_2nd));
	w_b.copyTo(w_pad_2nd(ranges_2nd));

	// Precalculating linear index steps for each new padded dimension
	const unsigned int step_y_pad_2nd = (width_ + padding * 2) * range_;
	const unsigned int step_x_pad_2nd = range_;
	const unsigned int step_z_pad_2nd = 1; // Unchanged

	if (DEBUG) {
		if (w_i_pad_2nd.isContinuous() && w_pad_2nd.isContinuous() &&
			w_i_b.isContinuous() && w_b.isContinuous()) {
			std::cout << "with continuous matrices" << std::endl;
		}
		else {
			std::cout << std::endl;
		}
	}

	for (int z = 0; z < range_; ++z) {
		for (int y = 0; y < height_; ++y) {
			for (int x = 0; x < width_; ++x) {
				Vec<double, KERNEL_SIZE> w_i_elems, w_elems;

				if (w_i_pad_2nd.isContinuous()) {
					int linear_idx_pad = y*step_y_pad_2nd + x*step_x_pad_2nd + z*step_z_pad_2nd;
					int linear_idx = y*step_y + x*step_x + z*step_z;

					w_i_elems = {
						w_i_pad_2nd_p[linear_idx_pad + 0 * step_x_pad_2nd],
						w_i_pad_2nd_p[linear_idx_pad + 1 * step_x_pad_2nd],
						w_i_pad_2nd_p[linear_idx_pad + 2 * step_x_pad_2nd],
						w_i_pad_2nd_p[linear_idx_pad + 3 * step_x_pad_2nd],
						w_i_pad_2nd_p[linear_idx_pad + 4 * step_x_pad_2nd]
					};

					w_i_b_p[linear_idx] = kernel.dot(w_i_elems);

					w_elems = {
						w_pad_2nd_p[linear_idx_pad + 0 * step_x_pad_2nd],
						w_pad_2nd_p[linear_idx_pad + 1 * step_x_pad_2nd],
						w_pad_2nd_p[linear_idx_pad + 2 * step_x_pad_2nd],
						w_pad_2nd_p[linear_idx_pad + 3 * step_x_pad_2nd],
						w_pad_2nd_p[linear_idx_pad + 4 * step_x_pad_2nd]
					};

					w_b_p[linear_idx] = kernel.dot(w_elems);

				}
				else {
					int idx1[3] = { y, x + 0, z };
					int idx2[3] = { y, x + 1, z };
					int idx3[3] = { y, x + 2, z };
					int idx4[3] = { y, x + 3, z };
					int idx5[3] = { y, x + 4, z };

					w_i_elems = {
						w_i_pad_2nd.at<double>(idx1),
						w_i_pad_2nd.at<double>(idx2),
						w_i_pad_2nd.at<double>(idx3),
						w_i_pad_2nd.at<double>(idx4),
						w_i_pad_2nd.at<double>(idx5)
					};

					w_i_b.at<double>(idx1) = kernel.dot(w_i_elems);

					w_elems = {
						w_pad_2nd.at<double>(idx1),
						w_pad_2nd.at<double>(idx2),
						w_pad_2nd.at<double>(idx3),
						w_pad_2nd.at<double>(idx4),
						w_pad_2nd.at<double>(idx5)
					};

					w_b.at<double>(idx1) = kernel.dot(w_elems);
				} // if - else
			} // 3rd dimension
		} // 2nd dimension
	} // 1st dimension

	  /*for (uchar i = 0; i < range_; ++i) {
		  Mat slice;
		  extractZSliceOBL(w_i_b.mul(.0001), i, slice);
		  imshow("w_i_b_2nd", slice);
		  waitKey(0);
	  }*/

	  // Deleting w_i_pad_2nd and w_pad_2nd to free memory
	w_i_pad_2nd.release();
	w_pad_2nd.release();

	// Convolving along the 1st dimension
	if (DEBUG) std::cout << "DEBUG: Convolving along the 1st dimension ";
	// Padding image for the convolution
	Vec<int, 3> padding_1st = { padding, 0, 0 };
	Range ranges_1st[3] = {
		Range(padding, padding + height_), Range::all(), Range::all()
	};
	Mat w_i_pad_1st = Mat(3, (int*)Mat(padding_1st + sizes + padding_1st).data, CV_64F, Scalar(0));
	Mat w_pad_1st = Mat(3, (int*)Mat(padding_1st + sizes + padding_1st).data, CV_64F, Scalar(0));
	double *w_i_pad_1st_p = (double*)w_i_pad_1st.data;
	double *w_pad_1st_p = (double*)w_pad_1st.data;

	w_i_b.copyTo(w_i_pad_1st(ranges_1st));
	w_b.copyTo(w_pad_1st(ranges_1st));

	// No need to calculate new linear indeces steps 
	// as padding the first dimension does not affect them at all

	if (DEBUG) {
		if (w_i_pad_1st.isContinuous() && w_pad_1st.isContinuous() &&
			w_i_b.isContinuous() && w_b.isContinuous()) {
			std::cout << "with continuous matrices" << std::endl;
		}
		else {
			std::cout << std::endl;
		}
	}
	for (int z = 0; z < range_; ++z) {
		for (int x = 0; x < width_; ++x) {
			for (int y = 0; y < height_; ++y) {
				Vec<double, KERNEL_SIZE> w_i_elems, w_elems;

				if (w_i_pad_1st.isContinuous()) {
					int linear_idx = y*step_y + x*step_x + z*step_z;

					w_i_elems = {
						w_i_pad_1st_p[linear_idx + 0 * step_y],
						w_i_pad_1st_p[linear_idx + 1 * step_y],
						w_i_pad_1st_p[linear_idx + 2 * step_y],
						w_i_pad_1st_p[linear_idx + 3 * step_y],
						w_i_pad_1st_p[linear_idx + 4 * step_y]
					};

					w_i_b_p[linear_idx] = kernel.dot(w_i_elems);

					w_elems = {
						w_pad_1st_p[linear_idx + 0 * step_y],
						w_pad_1st_p[linear_idx + 1 * step_y],
						w_pad_1st_p[linear_idx + 2 * step_y],
						w_pad_1st_p[linear_idx + 3 * step_y],
						w_pad_1st_p[linear_idx + 4 * step_y]
					};

					w_b_p[linear_idx] = kernel.dot(w_elems);

				}
				else {
					int idx1[3] = { y + 0, x, z };
					int idx2[3] = { y + 1, x, z };
					int idx3[3] = { y + 2, x, z };
					int idx4[3] = { y + 3, x, z };
					int idx5[3] = { y + 4, x, z };

					w_i_elems = {
						w_i_pad_1st.at<float>(idx1),
						w_i_pad_1st.at<float>(idx2),
						w_i_pad_1st.at<float>(idx3),
						w_i_pad_1st.at<float>(idx4),
						w_i_pad_1st.at<float>(idx5)
					};

					w_i_b.at<double>(idx1) = kernel.dot(w_i_elems);

					w_elems = {
						w_pad_1st.at<float>(idx1),
						w_pad_1st.at<float>(idx2),
						w_pad_1st.at<float>(idx3),
						w_pad_1st.at<float>(idx4),
						w_pad_1st.at<float>(idx5)
					};

					w_b.at<double>(idx1) = kernel.dot(w_elems);
				} // if - else
			} // 3rd dimension
		} // 2nd dimension
	} // 1st dimension

	  //for (uchar i = 0; i < range_; ++i) {
		 // Mat slice1, slice2;
		 // extractZSliceOBL(w_b.mul(.01), i, slice1);
		 // extractZSliceOBL(w_i_b.mul(.0001), i, slice2);
		 // //imshow("slice1", slice1);
		 // imshow("w_i_b", slice2);
		 // waitKey(0);
	  //}

	  // Deleting w_i_pad_1st and w_pad_1st to free memory
	w_i_pad_1st.release();
	w_pad_1st.release();

	Mat dst = _dst.getMat();
	uchar *dst_p = dst.data;

	// Interpolating to get the filtered image
	if (DEBUG) {
		if (w_i_b.isContinuous() && w_b.isContinuous() &&
			dst.isContinuous()) {
			std::cout << "DEBUG: Interpolating with continuous matrices" << std::endl;
		}
		else {
			std::cout << "DEBUG: Interpolating" << std::endl;
		}
	}

	// Padding arrays for interpolation
	Vec<int, 3> padding_all = { 1, 1, 1 };
	Range ranges[3] = {
		Range(0, height_), Range(0, width_), Range(0, range_)
	};

	Mat w_i_b_pad = Mat(3, (int*)Mat(sizes + padding_all).data, CV_64F, Scalar(0));
	Mat w_b_pad = Mat(3, (int*)Mat(sizes + padding_all).data, CV_64F, Scalar(0));

	w_i_b.copyTo(w_i_b_pad(ranges));
	w_b.copyTo(w_b_pad(ranges));

	/*for (uchar i = 0; i < range_+1; ++i) {
	Mat slice1, slice2;
	extractZSliceOBL(w_i_b_pad.mul(.0001), i, slice1);
	extractZSliceOBL(w_b_pad.mul(.01), i, slice2);
	imshow("slice1", slice1);
	imshow("slice2", slice2);
	waitKey(0);
	}*/

	for (uint Y = 0; Y < height; ++Y) {
		for (uint X = 0; X < width; ++X) {

			double qy = Y / sigmaSpatial;
			double qx = X / sigmaSpatial;
			double qz = (src.isContinuous() ?
				image_p[Y*width + X] : src.at<uchar>(Y, X)) / sigmaRange;

			double WI = interpolate(w_i_b_pad, qy, qx, qz);
			double W = interpolate(w_b_pad, qy, qx, qz);

			if (dst.isContinuous()) {
				dst_p[Y*width + X] = static_cast<uchar>(WI / W);
				//dst_p[Y*width + X] = (WI / W) > 255 ? 255 : (WI / W);
			}
			else {
				dst.at<uchar>(Y, X) = static_cast<uchar>(WI / W);
				//dst.at<uchar>(Y, X) = (WI / W) > 255 ? 255 : (WI / W);
			}
		}
	}

	t = ((double)getTickCount() - t) / getTickFrequency();

	if (DEBUG) {
		std::cout << "DEBUG: Image filtered in " << t << " seconds" << std::endl;
		//getchar();
	}

	return t;
}