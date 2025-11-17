
#include <sys/time.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "timer.h"	// 计时头文件

// vitis hls 相关库文件
#include "ap_fixed.h"	//ap_fixed<18,6,AP_TRN_ZERO,AP_SAT>		<W,I,Q,O,N>
#include "ap_int.h"	//ap_int<N> or ap_uint<N>, 1<=N<=1024
#include "hls_math.h"	//data_t s = hls::sinf(angle);
#include "hls_stream.h"

// 重要宏定义 通过宏定义来切换不同功能。
#define DEBUG			// 调试 宏定义

#ifdef DEBUG
	#define DEBUG_LOG(x) std::cout << "[DEBUG] " << x << std::endl
	#define DEBUG_INFO(x) std::cout << "[INFO] " << x << std::endl
	#define DEBUG_ERROR(x) std::cout << "[ERROR] " << x << std::endl
	#define DEBUG_TIME(x) std::cout << "[TIME] " << x << std::endl
	#define DEBUG_OPERATION(x) x
	#define DEBUG_GETCHAR
#endif

typedef float type_point;
typedef ap_fixed<20, 10> type_point_hw;
// typedef float type_point_hw;

typedef ap_fixed<20, 10> type_dist_hw;
// typedef ap_fixed<10, 4> type_dist_hw;
// typedef float type_dist_hw;

typedef ap_uint<32> voxel_int;
typedef ap_uint<8> flag_int;

struct My_PointXYZI {
    type_point x;
    type_point y;
    type_point z;
    // float intensity;
};

typedef struct PointInformation {
    // int hash;
    int encode;
    int flag;
    int depth;
} point_information;

static const int k_div_flag[4] = {
    3,          // k_div_flag[0] Dsvs
    7,         // k_div_flag[1] Oc1
    9,        // k_div_flag[2] Oc2
    11        // k_div_flag[3] Oc3
};
static const int k_no_div_flag = 27;

// 参数定义区
static const type_point_hw k_DSVS_split_unit_hw    = 2         ;
static const type_point_hw k_mid_split_unit_hw     = 1         ;
static const type_point k_DSVS_split_unit    = 2         ;
static const type_point k_mid_split_unit     = 1         ;

static const int K     = 1         ;
// static const type_dist_hw dist_unit_hw  = 0.2       ; //&&0.4不掉精度
static const int k_voxels_number_max    = 163400     ;
// 1-1238230; 2-163400; 3-51520; 3.4-33936; 3.6-30528; 4-20640; 4.2-18860; 4.4-17160; 4.6-12600; 5-10488
static const int k_max_encode           = 64        ;//512-4096

static const int k_Octree_threshold     = 64         ;// 找最优
static const int k_max_point_num        = 64         ; // 每个查询点的候选点最大数量4096 2048

static const int k_Depth                = 8         ;
static const int k_encode               = 8         ;
static const int k_reference_set_size 	= 1100000   ;// 1100000 - 3051163 - 7000000
static const int k_query_set_size 		= 10296     ;// 1001 - 3000 - 7296 - 10296
static const int k_data_max_value_abs 	= 500       ;
static const int k_axis_voxel_max 		= 500       ;

// 函数区
// 功能函数
void read_points_from_txt(std::string file_name, My_PointXYZI* laserCloudInArray, int & point_size);
void brute_force_search(My_PointXYZI* query_set, int query_set_size, My_PointXYZI* reference_set, int reference_set_size, My_PointXYZI* bf_KNN_query_result);
void compare_result(My_PointXYZI* query_set, int query_set_size, My_PointXYZI* new_query_result, My_PointXYZI* gt_query_result, int & error_count);

// 硬件初始化
void get_min_max(My_PointXYZI* data_set, int data_set_size, type_point_hw& data_set_xmin_hw, type_point_hw& data_set_ymin_hw, type_point_hw& data_set_zmin_hw);
void split_voxel_AABB(type_point split_unit, int& split_x_size_PL_hw, int& split_y_size_PL_hw, int& split_z_size_PL_hw, int& total_voxel_size_PL_hw);
void setup_hardware_PL(My_PointXYZI* KNN_reference_set, int reference_set_size,
    type_point_hw& data_set_xmin_hw, type_point_hw& data_set_ymin_hw, type_point_hw& data_set_zmin_hw, type_point_hw& voxel_split_unit_hw,
    int& split_x_size_PL_hw, int& split_y_size_PL_hw, int& split_z_size_PL_hw, int& total_voxel_size_PL_hw);

// DSVS Init
void DSVS_count_point(int* data_hash, int* DSVS_point_num, int KNN_reference_set_size);
void DSVS_cal_first_index(int* DSVS_point_num, int* DSVS_first_index, int total_voxel_size);
void DSVS_cal_flag(int* DSVS_first_index, int* DSVS_voxel_flag);
void DSVS_reorder(My_PointXYZI* set, int set_size, int* data_hash, int* DSVS_first_index, int* DSVS_point_num, int* count_point, 
    int* ordered_index, My_PointXYZI* ordered_point);
void DSVS_build(My_PointXYZI* KNN_reference_set, int reference_set_size, int* DSVS_point_num, int* DSVS_first_index, 
    int* DSVS_voxel_flag, My_PointXYZI* DSVS_ordered_set, int* DSVS_ordered_index);

// Octree Init
void dep1_Oc_encode(My_PointXYZI* ordered_point, int* DSVS_voxel_flag, int* DSVS_first_index, int *DSVS_point_num, point_information* point_inform);
void dep1_Oc_point_num(int* Oc_voxel_point_num, int* Oc_voxel_depth, int* DSVS_voxel_flag, int* DSVS_first_index, 
    int *DSVS_point_num, point_information* point_inform);
void cal_Oc_first_index(int Depth, int* Oc_voxel_flag, int* Oc_voxel_first_index, int* Oc_voxel_point_num);
void cal_Oc_voxel_flag(int Depth, int* Oc_voxel_flag, int* Oc_voxel_depth, int* Oc_voxel_point_num);
void cal_Oc_encode(int Depth, int* Oc_voxel_flag, int* Oc_voxel_first_index, int* Oc_voxel_point_num,
    My_PointXYZI* ordered_point, point_information* point_inform);
void depN_cal_Oc_voxel_point(int Depth, int* Oc_voxel_point_num, int* next_dep_point_num, int* Oc_voxel_first_index, 
    int* Oc_voxel_depth, int* Oc_voxel_flag, point_information* point_inform);
void depN_Octree_reorder(int Depth, My_PointXYZI* set, int* set_index, int set_size, int* data_hash, My_PointXYZI* depthN_set, int* depthN_set_index,
    int* Oc_voxel_first_index, point_information* point_inform, point_information* depN_point_inform, int* count_temp);
void Octree_Init(My_PointXYZI* DSVS_ordered_set, int* DSVS_ordered_index, int set_size, int* DSVS_point_num, int* DSVS_first_index,
    int* DSVS_voxel_flag, int* Oc_voxel_depth, int* Oc_voxel_point_num, int* Oc_voxel_flag, int* Oc_voxel_first_index,
    point_information* point_inform, point_information* dep1_point_inform, My_PointXYZI* dep1_set, int* dep1_set_index);
void Octree_build(int Depth, My_PointXYZI* dep1_set, int* dep1_set_index, int set_size, int* Oc_voxel_depth, int* Oc_voxel_point_num,
    int* Oc_voxel_flag, int* Oc_voxel_first_index, point_information* dep1_point_inform, point_information* dep2_point_inform,
    int* dep2_Oc_point_num, My_PointXYZI* dep2_set, int* dep2_set_index);
void Octree_test(My_PointXYZI* DSVS_ordered_set, int* DSVS_ordered_index, int set_size, int* DSVS_point_num, int* DSVS_first_index,
    int* DSVS_voxel_flag, int* Oc_voxel_depth, int* Oc_depN_point_num, int* Oc_voxel_flag, int* Oc_voxel_first_index, 
    point_information* point_inform, point_information* depN_point_inform, My_PointXYZI* depN_set, int* depN_set_index);

// Query
void cal_query_point_inform_new(My_PointXYZI* query_set, int q_set_size, int* q_hash, int* DSVS_voxel_flag, int* Oc_voxel_flag, 
    int* DSVS_first_index, int* Oc_voxel_first_index, int* q_flag, voxel_int* q_set_start_index);
void query_init(My_PointXYZI* query_set, int q_set_size, int* DSVS_voxel_flag, int* DSVS_first_index, 
    int* Oc_voxel_first_index, int* Oc_voxel_flag, int* q_flag, voxel_int* q_set_start_index);

void query_init_v2(My_PointXYZI* query_set, int q_set_size, int* DSVS_voxel_flag, int* DSVS_first_index, 
    int* q_flag, voxel_int* q_set_start_index);
void cal_query_info_v2(int q_set_size, int* q_hash, int* DSVS_voxel_flag,
    int* DSVS_first_index, int* q_flag, voxel_int* q_set_start_index);
// HW区
struct My_PointXYZI_HW {
    type_point_hw x;
    type_point_hw y;
    type_point_hw z;
};

// 函数
type_dist_hw cal_dist_hw(My_PointXYZI_HW data1, My_PointXYZI_HW data2);
void input_src_hw(My_PointXYZI* set, hls::stream<My_PointXYZI_HW>& set_input, 
    voxel_int* q_set_start_index, hls::stream<voxel_int>& q_set_start_index_stream, 
    int* q_flag, hls::stream<int>& q_flag_stream, int set_size);
void compareSwap(type_dist_hw &a, type_dist_hw &b, voxel_int &a_index, voxel_int &b_index, bool dir);
void bitonicSort64_hw(type_dist_hw* arr, voxel_int* arr_index);
void bitonicSort128_hw(type_dist_hw* arr, voxel_int* arr_index);
void bitonicSort256_hw(type_dist_hw* arr, voxel_int* arr_index);
void bitonicSort512_hw(type_dist_hw *arr, voxel_int *arr_index);
void bitonicSort_hw(type_dist_hw* arr, voxel_int* arr_index, int num);
void quick_sort_hw(type_dist_hw* arr, voxel_int* arr_index);

void calculate_hash_hw(hls::stream<My_PointXYZI_HW>& set_input, hls::stream<My_PointXYZI_HW>& set_output, hls::stream<voxel_int>& data_set_hash, int set_size);
void cal_position_hw(hls::stream<My_PointXYZI_HW>& set_input, hls::stream<My_PointXYZI_HW>& set_output, hls::stream<voxel_int>& q_hash,
    hls::stream<int>& q_flag_stream, hls::stream<voxel_int>& q_start_index, hls::stream<voxel_int>& q_start_index_out, 
    int* Oc_voxel_flag, int* Oc_voxel_first_index, int q_set_size);
void cal_candidate_point_hw(hls::stream<My_PointXYZI_HW>& candidate_point_stream, hls::stream<voxel_int>& q_set_start_index, 
    hls::stream<voxel_int>& q_hash, hls::stream<bool>& update_flag_stream, My_PointXYZI_HW* Oc_ordered_set, int q_set_size);
void KNN_search_hw(My_PointXYZI_HW* result, hls::stream<My_PointXYZI_HW>& candidate_point_stream, 
    hls::stream<bool>& update_flag_stream, hls::stream<My_PointXYZI_HW>& q_set, int q_set_size);
void Ds_Oc_KNN_search_hw(My_PointXYZI_HW* result_hw, My_PointXYZI* query_set, voxel_int* q_set_start_index, int* q_flag,
    My_PointXYZI_HW* Oc_ordered_set, int* Oc_voxel_flag, int* Oc_voxel_first_index, int q_set_size,
    type_point_hw data_set_xmin_PL_hw, type_point_hw data_set_ymin_PL_hw, type_point_hw data_set_zmin_PL_hw, type_point_hw voxel_split_unit_PL_hw,
    voxel_int split_x_size_PL_hw, voxel_int split_y_size_PL_hw, voxel_int split_z_size_PL_hw, voxel_int total_voxel_size_PL_hw);

// bitonic512
void swap(type_dist_hw &a, type_dist_hw &b, voxel_int &a_idx, voxel_int &b_idx);
void compare_swap(type_dist_hw &a, type_dist_hw &b, voxel_int &a_idx, voxel_int &b_idx, bool dir);
void bitonic_merge(type_dist_hw arr[512], voxel_int arr_index[512], int low, int cnt, bool dir);
void bitonic_sort_rec(type_dist_hw arr[512], voxel_int arr_index[512], int low, int cnt, bool dir);
void bitonicSort512_hw(type_dist_hw arr[512], voxel_int arr_index[512]);
