#include <stdio.h>
#include <math.h>
#include <cstdlib>// Header file needed to use rand
#include <ctime> // Header file needed to use time
#include <iomanip>
// XRT includes
#include "xrt/xrt_bo.h"
#include <experimental/xrt_xclbin.h>
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"

#include "DSVS_Octree.h"

int main(int argc, char** argv)
{
    /****************************************** 0. 变量定义 ******************************/
    // ====================== 分配设备内存 ======================
    std::cout << "初始化 FPGA 设备, done! \n";
    int device_index = 0;
    std::string xclbin_path = "/home/wyt/workspack/vitis/DsOc_v2/hlsKernel/output/DsOc_kernel.xclbin";

    auto device = xrt::device(device_index);
    auto uuid = device.load_xclbin(xclbin_path);
    auto krnl = xrt::kernel(device, uuid, "Ds_Oc_KNN_search_hw");

    // 输入点云 参考点集 及 搜索集
    My_PointXYZI KNN_reference_set[k_reference_set_size];
    My_PointXYZI KNN_query_set[k_query_set_size];

    int reference_set_num = 0;          // 参考数据集帧索引，可通过命令行设置，默认0
    int query_set_num = 0;              // 搜索数据集帧索引，可通过命令行设置，默认0
    int reference_point_num = 1000;     // 参考点数量，单位k
    int query_point_num = 10;            // 查询点数量，单位k

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
    My_PointXYZI depN_set[KNN_reference_set_size];
    int depN_set_index[KNN_reference_set_size];

//query 变量
    int q_DSVS_point_num[k_voxels_number_max] = {0};
    int q_DSVS_first_index[k_voxels_number_max] = {0};
    int q_DSVS_voxel_flag[k_voxels_number_max] = {0};
    My_PointXYZI q_ordered_set[KNN_query_set_size];
    int q_ordered_index[KNN_query_set_size];

    voxel_int q_set_start_index[k_query_set_size];
    int q_flag[k_query_set_size];

// Parameter Init
    setup_hardware_PL(KNN_reference_set, KNN_reference_set_size, data_set_xmin_PL_hw, data_set_ymin_PL_hw, data_set_zmin_PL_hw, voxel_split_unit_PL_hw, split_x_size_PL_hw, split_y_size_PL_hw, split_z_size_PL_hw, total_voxel_size_PL_hw);
// DSVS Init
    DSVS_build(KNN_reference_set, KNN_reference_set_size, DSVS_point_num, DSVS_first_index, DSVS_voxel_flag, DSVS_ordered_set, DSVS_ordered_index);
// Octree Init
    Octree_test(DSVS_ordered_set, DSVS_ordered_index, KNN_reference_set_size, DSVS_point_num, DSVS_first_index, DSVS_voxel_flag, 
        Oc_voxel_depth, Oc_depN_point_num, Oc_voxel_flag, Oc_voxel_first_index, point_inform, depN_point_inform, depN_set, depN_set_index);
// Query Init
    DSVS_build(KNN_query_set, KNN_query_set_size, q_DSVS_point_num, q_DSVS_first_index, q_DSVS_voxel_flag, q_ordered_set, q_ordered_index);
    // query_init(q_ordered_set, KNN_query_set_size, DSVS_voxel_flag, DSVS_first_index, Oc_voxel_first_index, Oc_voxel_flag, q_flag, q_set_start_index);
    query_init_v2(q_ordered_set, KNN_query_set_size, DSVS_voxel_flag, DSVS_first_index, q_flag, q_set_start_index);
TIMER_STOP_ID(0);
DEBUG_TIME("Finish software building DSVS_Octree with reference points " << KNN_reference_set_size << " with " << TIMER_REPORT_MS(0) << " ms !" );
    
    DEBUG_LOG("           k_max_encode          : " << k_max_encode             );
    DEBUG_LOG("         k_Octree_threshold      : " << k_Octree_threshold       );
    DEBUG_LOG("        k_DSVS_split_unit_hw     : " << k_DSVS_split_unit_hw     );
    DEBUG_LOG("         split_x_size_PL_hw      : " << split_x_size_PL_hw       );
    DEBUG_LOG("         split_y_size_PL_hw      : " << split_y_size_PL_hw       );
    DEBUG_LOG("         split_z_size_PL_hw      : " << split_z_size_PL_hw       );
    DEBUG_LOG("  hash = total_voxel_size_PL_hw  : " << total_voxel_size_PL_hw   );

    My_PointXYZI_HW depN_set_hw[k_reference_set_size];
    for(int i = 0; i < KNN_reference_set_size; i++)
    {
        My_PointXYZI_HW temp;
        temp.x = (type_point_hw)depN_set[i].x;
        temp.y = (type_point_hw)depN_set[i].y;
        temp.z = (type_point_hw)depN_set[i].z;
        depN_set_hw[i] = temp;
    }
    /********************************************** 3. KNN 搜索 *******************************************************/

    std::cout << "Allocate Buffer in Global Memory\n";
    // size_t q_size_My_Point_HW   = sizeof(My_PointXYZI_HW) * k_query_set_size * K;
    size_t q_size_My_Point_HW   = sizeof(My_PointXYZI_HW) * k_query_set_size;
    size_t q_size_My_Point      = sizeof(My_PointXYZI) * k_query_set_size;
    size_t q_size_start_index   = sizeof(voxel_int) * k_query_set_size;
    size_t q_size_voxel_flag    = sizeof(int) * k_query_set_size;
    size_t ref_size_My_Point_HW = sizeof(My_PointXYZI_HW) * k_reference_set_size;
    size_t VOS_info_size        = sizeof(int) * k_voxels_number_max * 64;

    auto bo_result_hw           = xrt::bo(device, q_size_My_Point_HW, krnl.group_id(0));
    auto bo_query_set           = xrt::bo(device, q_size_My_Point, krnl.group_id(1));
    auto bo_start_index         = xrt::bo(device, q_size_start_index, krnl.group_id(2));
    auto bo_q_flag              = xrt::bo(device, q_size_voxel_flag, krnl.group_id(3));
    auto bo_ref_set             = xrt::bo(device, ref_size_My_Point_HW, krnl.group_id(4));
    auto bo_VOS_flag            = xrt::bo(device, VOS_info_size, krnl.group_id(5));
    auto bo_VOS_index           = xrt::bo(device, VOS_info_size, krnl.group_id(6));

    // ====================== 映射主机内存 ======================
    auto bo0_map                = bo_result_hw.map<My_PointXYZI_HW*>();
    auto bo1_map                = bo_query_set.map<My_PointXYZI*>();
    auto bo2_map                = bo_start_index.map<voxel_int*>();
    auto bo3_map                = bo_q_flag.map<int*>();
    auto bo4_map                = bo_ref_set.map<My_PointXYZI_HW*>();
    auto bo5_map                = bo_VOS_flag.map<int*>();
    auto bo6_map                = bo_VOS_index.map<int*>();

    // std::fill(bo0_map, bo0_map + k_query_set_size * K, My_PointXYZI_HW{0, 0, 0});
    std::fill(bo0_map, bo0_map + k_query_set_size, My_PointXYZI_HW{0, 0, 0});
    std::fill(bo1_map, bo1_map + k_query_set_size, My_PointXYZI{0, 0, 0});
    std::fill(bo2_map, bo2_map + k_query_set_size, 0);
    std::fill(bo3_map, bo3_map + k_query_set_size, 0);
    std::fill(bo4_map, bo4_map + k_reference_set_size, My_PointXYZI_HW{0, 0, 0});
    std::fill(bo5_map, bo5_map + k_voxels_number_max * 64, 0);
    std::fill(bo6_map, bo6_map + k_voxels_number_max * 64, 0);

    std::cout << "Finish Allocate Buffer in Global Memory\n";

    for(int i = 0; i < k_query_set_size; i++)
    {
        bo1_map[i] = KNN_query_set[i];
        bo2_map[i] = q_set_start_index[i];
        bo3_map[i] = q_flag[i];
    }
    for(int i = 0; i < k_reference_set_size; i++)
    {
        bo4_map[i] = depN_set_hw[i];
    }
    for(int i = 0; i < k_voxels_number_max * 64; i++)
    {
        bo5_map[i] = Oc_voxel_flag[i];
        bo6_map[i] = Oc_voxel_first_index[i];
    }
    
    bo_query_set.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_start_index.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_q_flag.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_ref_set.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_VOS_flag.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    bo_VOS_index.sync(XCL_BO_SYNC_BO_TO_DEVICE);

TIMER_START(2);
    // 记录开始时间点
    auto start_time = std::chrono::high_resolution_clock::now();
    // 调用内核函数
    auto run = krnl(bo_result_hw, bo_query_set, bo_start_index, bo_q_flag, bo_ref_set, bo_VOS_flag, bo_VOS_index, KNN_query_set_size,
        data_set_xmin_PL_hw, data_set_ymin_PL_hw, data_set_zmin_PL_hw, voxel_split_unit_PL_hw,
        split_x_size_PL_hw, split_y_size_PL_hw, split_z_size_PL_hw, total_voxel_size_PL_hw);
    run.wait();
    // 计算耗时
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_ms = end_time - start_time;
    std::cout << "Kernel execution time: " << std::fixed << std::setprecision(3) << duration_ms.count() << " ms" << std::endl;

    
TIMER_STOP_ID(2);
DEBUG_TIME("Finish search DSVS with query points " << KNN_query_set_size << " with " << TIMER_REPORT_MS(2) << " ms !" );

    bo_result_hw.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    // 将结果从设备内存复制到主机内存
    // My_PointXYZI_HW result_hw[k_query_set_size * K];
    // for(int i = 0; i < k_query_set_size * K; i++)
    My_PointXYZI_HW result_hw[k_query_set_size];
    for(int i = 0; i < k_query_set_size; i++)
    {
        result_hw[i] = bo0_map[i];
    }


// // 暴力搜索~~~~！！！！我突然发现我好像不能这么比，因为我是类型转换的，所以下面也得用类型转换的数据比较
//     My_PointXYZI bf_KNN_query_result[KNN_query_set_size];
//     brute_force_search(KNN_query_set, KNN_query_set_size, KNN_reference_set, KNN_reference_set_size, bf_KNN_query_result);

//     My_PointXYZI result[k_query_set_size];
//     for(int i = 0; i < k_query_set_size; i++)
//     {
//         My_PointXYZI temp;
//         temp.x = (float)result_hw[i].x;
//         temp.y = (float)result_hw[i].y;
//         temp.z = (float)result_hw[i].z;
//         result[i] = temp;
//     }

//     /********************************************** 4. KNN 结果比较 *******************************************************/
//     int error_count = 0;
//     float NN_true_ratio;
// // 比较结果
//     compare_result(KNN_query_set, KNN_query_set_size, result, bf_KNN_query_result, error_count);
//     DEBUG_LOG("image_KNN VS brute-force: " << ", error_count: " << error_count << " with true ratio = " << float(1-(float)error_count/(float)KNN_query_set_size));

//     NN_true_ratio = float(1-(float)error_count/(float)KNN_query_set_size);

//      /********************************************** 5. 验证结果 *******************************************************/
//      if (NN_true_ratio < 0.01)
//      {
//          DEBUG_INFO("NN_true_ratio : " << NN_true_ratio);
//          fprintf(stdout, "*******************************************\n");
//          fprintf(stdout, "FAIL: Output DOES NOT match the golden output\n");
//          fprintf(stdout, "*******************************************\n");
//          return 1;
//      } 
//      else 
//      {
//          DEBUG_INFO("NN_true_ratio: " << NN_true_ratio);
//          fprintf(stdout, "*******************************************\n");
//          fprintf(stdout, "PASS: The output matches the golden output!\n");
//          fprintf(stdout, "*******************************************\n");
//          return 0;
//      }

}
