#include <stdio.h>
#include <math.h>
#include <cstdlib>// Header file needed to use rand
#include <ctime> // Header file needed to use time
#include <iomanip>

#include "DSVS_Octree.h"

int main(int argc, char** argv)
{
    /****************************************** 0. 变量定义 ******************************/
    // 输入点云 参考点集 及 搜索集
    My_PointXYZI KNN_reference_set[k_reference_set_size];
    My_PointXYZI KNN_query_set[k_query_set_size];

    int reference_set_num = 0;          // 参考数据集帧索引，可通过命令行设置，默认0
    int query_set_num = 0;              // 搜索数据集帧索引，可通过命令行设置，默认0
    int reference_point_num = 1000;     // 参考点数量，单位k
    int query_point_num = 1;            // 查询点数量，单位k

    std::string pointcloud_dir = "/home/wyt/Documents/Data/fpga_map_knn_test_data/";

    std::cout << "Avoid Segmentation fault (core dumped) : ulimit -s unlimited" << std::endl;

    DEBUG_LOG("Input pointcloud dir: " << pointcloud_dir);
    
    /****************************************** 1. 从文件中读取点云数据 ******************************/
    std::string KNN_reference_set_file_name;
    std::string KNN_query_set_file_name;
    
    KNN_reference_set_file_name = pointcloud_dir + "dataset_map_" + std::to_string(reference_set_num) + "_" + std::to_string(reference_point_num) + "k_points.txt";
    KNN_query_set_file_name = pointcloud_dir + "dataset_query_" + std::to_string(query_set_num)  + "_" + std::to_string(query_point_num) + "k_points.txt";

    DEBUG_INFO("KNN_reference_set_file_name: " << KNN_reference_set_file_name);
    DEBUG_INFO("KNN_query_set_file_name: " << KNN_query_set_file_name);

    // 从点云文件读取点云
    int KNN_reference_set_size, KNN_query_set_size;
    DEBUG_LOG("reading KNN_reference_set ......");
    read_points_from_txt(KNN_reference_set_file_name, KNN_reference_set, KNN_reference_set_size);
    read_points_from_txt(KNN_query_set_file_name, KNN_query_set, KNN_query_set_size);

    DEBUG_LOG("KNN_query_set_size: " << KNN_query_set_size);
    DEBUG_LOG("KNN_reference_set_size: " << KNN_reference_set_size);

    /********************************************** 2. 建立 数据结构 *******************************************************/
TIMER_INIT(7);
TIMER_START(0);
// DSVS 变量
    type_point_hw data_set_xmin_PL_hw;
    type_point_hw data_set_ymin_PL_hw;
    type_point_hw data_set_zmin_PL_hw;
    type_point_hw voxel_split_unit_PL_hw;
    int split_x_size_PL_hw;
    int split_y_size_PL_hw;
    int split_z_size_PL_hw;
    int total_voxel_size_PL_hw;

    int DSVS_point_num[k_voxels_number_max] = {0};
    int DSVS_first_index[k_voxels_number_max] = {0};
    int DSVS_voxel_flag[k_voxels_number_max] = {0};
    My_PointXYZI DSVS_ordered_set[KNN_reference_set_size];
    int DSVS_ordered_index[KNN_reference_set_size];

// Octree 变量
    int Oc_voxel_depth[k_voxels_number_max * k_max_encode] = {0};//depth_int Oc_voxel_depth[k_voxels_number_max * k_max_encode] = {0};
    int Oc_voxel_flag[k_voxels_number_max * k_max_encode] = {0};
    int Oc_voxel_first_index[k_voxels_number_max * k_max_encode] = {0};
    int Oc_depN_point_num[k_voxels_number_max * k_max_encode] = {0};
    point_information point_inform[KNN_reference_set_size];
    point_information depN_point_inform[KNN_reference_set_size];
    My_PointXYZI dep2_set[KNN_reference_set_size];
    int dep2_set_index[KNN_reference_set_size];

// // Parameter Init
//     setup_hardware_PL(KNN_reference_set, KNN_reference_set_size, data_set_xmin_PL_hw, data_set_ymin_PL_hw, data_set_zmin_PL_hw, voxel_split_unit_PL_hw, split_x_size_PL_hw, split_y_size_PL_hw, split_z_size_PL_hw, total_voxel_size_PL_hw);
// // DSVS Init
//     DSVS_build(KNN_reference_set, KNN_reference_set_size, DSVS_point_num, DSVS_first_index, DSVS_voxel_flag, DSVS_ordered_set, DSVS_ordered_index);
// // Octree Init
//     Octree_test(DSVS_ordered_set, DSVS_ordered_index, KNN_reference_set_size, DSVS_point_num, DSVS_first_index, DSVS_voxel_flag, 
//         Oc_voxel_depth, Oc_depN_point_num, Oc_voxel_flag, Oc_voxel_first_index, point_inform, depN_point_inform, dep2_set, dep2_set_index);
    
// TIMER_STOP_ID(0);
// DEBUG_TIME("Finish software building DSVS_Octree with reference points " << KNN_reference_set_size << " with " << TIMER_REPORT_MS(0) << " ms !" );
    
//     DEBUG_LOG("           k_max_encode          : " << k_max_encode             );
//     DEBUG_LOG("         k_Octree_threshold      : " << k_Octree_threshold       );
//     DEBUG_LOG("        k_DSVS_split_unit_hw     : " << k_DSVS_split_unit_hw     );
//     DEBUG_LOG("         split_x_size_PL_hw      : " << split_x_size_PL_hw       );
//     DEBUG_LOG("         split_y_size_PL_hw      : " << split_y_size_PL_hw       );
//     DEBUG_LOG("         split_z_size_PL_hw      : " << split_z_size_PL_hw       );
//     DEBUG_LOG("  hash = total_voxel_size_PL_hw  : " << total_voxel_size_PL_hw   );

//     /********************************************** 3. KNN 搜索 *******************************************************/

    type_dist_hw distance[k_max_point_num] = {
    173, 26, 31, 130, 191, 45, 168, 25, 79, 55, 153, 249, 70, 43, 9, 235,
    214, 156, 193, 186, 37, 132, 205, 1, 66, 208, 93, 247, 4, 239, 145, 252,
    94, 199, 29, 118, 78, 182, 183, 211, 12, 167, 144, 184, 115, 61, 111, 83,
    67, 240, 200, 217, 96, 158, 87, 244, 0, 135, 206, 14, 254, 123, 51, 184,
    41, 35, 44, 169, 215, 237, 232, 3, 64, 113, 6, 124, 221, 21, 172, 28,
    112, 141, 46, 92, 163, 138, 131, 198, 7, 80, 184, 34, 53, 131, 104, 97,
    38, 76, 2, 36, 90, 177, 48, 159, 153, 152, 216, 65, 72, 10, 42, 147,
    110, 15, 180, 195, 136, 23, 224, 85, 113, 88, 227, 122, 107, 99, 33, 158,
    118, 162, 134, 56, 16, 71, 27, 246, 130, 121, 175, 74, 66, 64, 225, 102,
    174, 114, 119, 254, 207, 55, 91, 100, 181, 5, 95, 109, 84, 177, 222, 187,
    190, 40, 242, 103, 30, 223, 158, 20, 239, 198, 60, 157, 139, 68, 150, 54,
    197, 146, 19, 177, 98, 69, 77, 196, 178, 11, 179, 43, 170, 220, 126, 134,
    132, 118, 167, 188, 73, 25, 125, 13, 117, 162, 136, 17, 232, 181, 52, 114,
    184, 75, 22, 189, 39, 201, 181, 24, 50, 150, 226, 229, 105, 108, 148, 49,
    137, 234, 223, 206, 248, 217, 138, 105, 176, 165, 218, 117, 43, 193, 202, 58,
    254, 128, 162, 42, 9, 106, 64, 53, 154, 151, 157, 220, 62, 140, 237, 41
};
    voxel_int index[k_max_point_num] = {
    173, 26, 31, 130, 191, 45, 168, 25, 79, 55, 153, 249, 70, 43, 9, 235,
    214, 156, 193, 186, 37, 132, 205, 1, 66, 208, 93, 247, 4, 239, 145, 252,
    94, 199, 29, 118, 78, 182, 183, 211, 12, 167, 144, 184, 115, 61, 111, 83,
    67, 240, 200, 217, 96, 158, 87, 244, 0, 135, 206, 14, 254, 123, 51, 184,
    41, 35, 44, 169, 215, 237, 232, 3, 64, 113, 6, 124, 221, 21, 172, 28,
    112, 141, 46, 92, 163, 138, 131, 198, 7, 80, 184, 34, 53, 131, 104, 97,
    38, 76, 2, 36, 90, 177, 48, 159, 153, 152, 216, 65, 72, 10, 42, 147,
    110, 15, 180, 195, 136, 23, 224, 85, 113, 88, 227, 122, 107, 99, 33, 158,
    118, 162, 134, 56, 16, 71, 27, 246, 130, 121, 175, 74, 66, 64, 225, 102,
    174, 114, 119, 254, 207, 55, 91, 100, 181, 5, 95, 109, 84, 177, 222, 187,
    190, 40, 242, 103, 30, 223, 158, 20, 239, 198, 60, 157, 139, 68, 150, 54,
    197, 146, 19, 177, 98, 69, 77, 196, 178, 11, 179, 43, 170, 220, 126, 134,
    132, 118, 167, 188, 73, 25, 125, 13, 117, 162, 136, 17, 232, 181, 52, 114,
    184, 75, 22, 189, 39, 201, 181, 24, 50, 150, 226, 229, 105, 108, 148, 49,
    137, 234, 223, 206, 248, 217, 138, 105, 176, 165, 218, 117, 43, 193, 202, 58,
    254, 128, 162, 42, 9, 106, 64, 53, 154, 151, 157, 220, 62, 140, 237, 41
};
    
    bitonicSort256_hw(distance, index);
    DEBUG_LOG(distance[0]);


}