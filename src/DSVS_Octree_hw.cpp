#include "DSVS_Octree.h"

static voxel_int split_x_array_size_hw;
static voxel_int split_y_array_size_hw;
static voxel_int split_z_array_size_hw;
static voxel_int total_voxel_size_hw;

static type_point_hw data_set_xmin_hw;
static type_point_hw data_set_ymin_hw;
static type_point_hw data_set_zmin_hw;
static type_point_hw voxel_split_unit_hw;

type_dist_hw cal_dist_hw(My_PointXYZI_HW data1, My_PointXYZI_HW data2) //data1是查询点，data2是参考点
{
#pragma HLS INLINE
    return ((data1.x - data2.x) * (data1.x - data2.x)
    + (data1.y - data2.y) * (data1.y - data2.y)
    + (data1.z - data2.z) * (data1.z - data2.z)
    );
}

void input_src_hw(My_PointXYZI* set, hls::stream<My_PointXYZI_HW>& set_input, 
    voxel_int* q_set_start_index, hls::stream<voxel_int>& q_set_start_index_stream, 
    int* q_flag, hls::stream<int>& q_flag_stream, int set_size)
{
    // DEBUG_OPERATION(int write_ref = 0;)
    for(int i = 0; i < set_size; i++)
    {
#pragma HLS pipeline II=1

        My_PointXYZI_HW temp_point_hw;
        voxel_int temp_index = q_set_start_index[i];
        int temp_flag = q_flag[i];

        temp_point_hw.x = set[i].x;
        temp_point_hw.y = set[i].y;
        temp_point_hw.z = set[i].z;
        set_input.write(temp_point_hw);
        q_set_start_index_stream.write(temp_index);
        q_flag_stream.write(temp_flag);
        // DEBUG_OPERATION(write_ref ++;)
    }
    // DEBUG_OPERATION(DEBUG_LOG("input_src_hw write_ref: " << write_ref);)
}

void calculate_hash_hw(hls::stream<My_PointXYZI_HW>& set_input, hls::stream<My_PointXYZI_HW>& set_output, hls::stream<voxel_int>& data_set_hash, int set_size)
{
    for (int i = 0; i < set_size; i++) 
    {
#pragma HLS LOOP_TRIPCOUNT max=k_query_set_size min=k_query_set_size
#pragma HLS pipeline II=1
        My_PointXYZI_HW temp_point = set_input.read();
        set_output.write(temp_point);

        type_point_hw data_x = temp_point.x;
        type_point_hw data_y = temp_point.y;
        type_point_hw data_z = temp_point.z;
        voxel_int x_split_array_size = split_x_array_size_hw;
        voxel_int y_split_array_size = split_y_array_size_hw;
        voxel_int z_split_array_size = split_z_array_size_hw;

        voxel_int x_index;
        voxel_int y_index;
        voxel_int z_index;

        if (data_x <= data_set_xmin_hw)
            x_index = 0;
        else
            x_index = (voxel_int)((data_x - data_set_xmin_hw) / voxel_split_unit_hw);

        if (x_index >= x_split_array_size)
            x_index = x_split_array_size - 1;

        if (data_y <= data_set_ymin_hw)
            y_index = 0;
        else
            y_index = (voxel_int)((data_y - data_set_ymin_hw) / voxel_split_unit_hw);

        if (y_index >= y_split_array_size)
            y_index = y_split_array_size - 1;

        if (data_z <= data_set_zmin_hw)
            z_index = 0;
        else
            z_index = (voxel_int)((data_z - data_set_zmin_hw) / voxel_split_unit_hw);

        if (z_index >= z_split_array_size)
            z_index = z_split_array_size - 1;

        //transform 3d index to a 1d index
        voxel_int data_hash = x_index * y_split_array_size * z_split_array_size + y_index * z_split_array_size + z_index;

        if (data_hash >= total_voxel_size_hw)
            data_hash = (total_voxel_size_hw - 1);
        if (data_hash < 0)
            data_hash = 0;
        
        data_set_hash.write(data_hash);
        // DEBUG_LOG("data_set_hash[" << i << "]: " << data_hash);
    }
}

// /*
void cal_position_hw(hls::stream<My_PointXYZI_HW>& set_input, hls::stream<My_PointXYZI_HW>& set_output, hls::stream<voxel_int>& q_hash,
    hls::stream<int>& q_flag_stream, hls::stream<voxel_int>& q_start_index, hls::stream<voxel_int>& q_start_index_out, 
    int* Oc_voxel_flag, int* Oc_voxel_first_index, int q_set_size)
{
    type_point_hw child_voxel_split_size = voxel_split_unit_hw / 2;
    for(int i = 0; i < q_set_size; i++)
    {
#pragma HLS LOOP_TRIPCOUNT max=k_query_set_size min=k_query_set_size
#pragma HLS pipeline II=1
        voxel_int temp_hash = q_hash.read();
        My_PointXYZI_HW temp_point = set_input.read();
        int temp_flag = q_flag_stream.read();
        voxel_int temp_index = q_start_index.read();
        set_output.write(temp_point);

        int x_gain = split_y_array_size_hw * split_z_array_size_hw;
        int x_index = temp_hash / x_gain;
        int y_gain = split_z_array_size_hw;
        int y_index = (temp_hash - x_index * x_gain) / y_gain;
        int z_index = temp_hash - x_index * x_gain - y_index * y_gain;

        type_point_hw xmin = data_set_xmin_hw + x_index * voxel_split_unit_hw;
        type_point_hw ymin = data_set_ymin_hw + y_index * voxel_split_unit_hw;
        type_point_hw zmin = data_set_zmin_hw + z_index * voxel_split_unit_hw;

        int code = 0;
        code |= (temp_point.x - xmin >= k_mid_split_unit_hw) << 2;
        code |= (temp_point.y - ymin >= k_mid_split_unit_hw) << 1;
        code |= (temp_point.z - zmin >= k_mid_split_unit_hw) << 0;

        if(temp_flag == k_div_flag[0])
        {
            temp_flag = Oc_voxel_flag[temp_hash*k_max_encode + code];
            temp_index = Oc_voxel_first_index[temp_hash*k_max_encode + code];
            if(temp_flag == k_div_flag[1])
            {
                type_point_hw xmin_child = xmin + (code & 0x4 ? child_voxel_split_size : type_point_hw(0));
                type_point_hw ymin_child = ymin + (code & 0x2 ? child_voxel_split_size : type_point_hw(0));
                type_point_hw zmin_child = zmin + (code & 0x1 ? child_voxel_split_size : type_point_hw(0));

                type_point_hw x_rel = (temp_point.x - xmin_child) / child_voxel_split_size;
                type_point_hw y_rel = (temp_point.y - ymin_child) / child_voxel_split_size;
                type_point_hw z_rel = (temp_point.z - zmin_child) / child_voxel_split_size;

                int sub_code = 0;
                sub_code |= (x_rel >= k_mid_split_unit_hw) << 2;
                sub_code |= (y_rel >= k_mid_split_unit_hw) << 1;
                sub_code |= (z_rel >= k_mid_split_unit_hw) << 0;

                int final_code = (code << 3) | sub_code;
                temp_index = Oc_voxel_first_index[temp_hash*k_max_encode + final_code];

                // if(sub_code == 0 | sub_code == 1 | sub_code == 2 | sub_code == 4)
                // {
                //     temp_index = temp_index - 12;
                // }
                // else
                // {
                //     temp_index = temp_index + 12;
                // }
            }
            // else
            // {
            //     if(code == 0 | code == 1 | code == 2 | code == 4)
            //     {
            //         // temp_index = temp_index - 12;
            //     }
            //     // else
            //     // {
            //     //     temp_index = temp_index + 12;
            //     // }
            // }
        }
        // if(code == 0 | code == 1 | code == 2 | code == 4)
        // {
        //     temp_index = temp_index - 12;
        // }
        // else
        // {
        //     temp_index = temp_index + 12;
        // }
        q_start_index_out.write(temp_index);
    }
}
// */

void cal_candidate_point_hw(hls::stream<My_PointXYZI_HW>& candidate_point_stream, hls::stream<voxel_int>& q_set_start_index, 
    hls::stream<bool>& update_flag_stream, My_PointXYZI_HW* Oc_ordered_set, int q_set_size)
{
    voxel_int previous_start_index = 0; // 用于存储上一个开始索引

    for(int i = 0; i < q_set_size; i++)
    {
#pragma HLS LOOP_TRIPCOUNT max=k_query_set_size min=k_query_set_size
#pragma HLS pipeline II=1

        voxel_int start_index = q_set_start_index.read();

        bool update_flag = false; // 新增标志位
        // 检查当前的 start_index 是否与上一个相同
        if (start_index != previous_start_index)
        {
            // 如果不同，继续正常选取候选点
            for(int j = 0; j < k_Octree_threshold; j++)
            {
// #pragma HLS UNROLL
#pragma HLS PIPELINE II=1
                candidate_point_stream.write(Oc_ordered_set[start_index + j]);
            }
            update_flag = true; // 设置更新标志位为真
        }
        else
        {
            // 如果相同，直接跳过选取过程
        }

        update_flag_stream.write(update_flag); // 将标志位写入流
        previous_start_index = start_index; // 更新之前的开始索引
    }
}

void KNN_search_hw(My_PointXYZI_HW* result, hls::stream<My_PointXYZI_HW>& candidate_point_stream, 
    hls::stream<bool>& update_flag_stream, hls::stream<My_PointXYZI_HW>& q_set, int q_set_size)
{
    type_dist_hw distance[k_max_point_num];
    voxel_int index[k_max_point_num];
    My_PointXYZI_HW candidate_buffer[k_max_point_num];

#pragma HLS bind_storage variable=distance type=ram_2p impl=bram
#pragma HLS bind_storage variable=index type=ram_2p impl=bram
// #pragma HLS bind_storage variable=candidate_buffer type=ram_2p impl=bram
// #pragma HLS ARRAY_PARTITION variable=distance complete
// #pragma HLS ARRAY_PARTITION variable=index complete

    for(int i = 0; i < q_set_size; i++)
    {
#pragma HLS PIPELINE II=1
        My_PointXYZI_HW q_temp_point = q_set.read();
        bool update_flag = update_flag_stream.read();

        if(update_flag)
        {
            for(int j = 0; j < k_Octree_threshold; j++)
            {
#pragma HLS PIPELINE II=1
// #pragma HLS UNROLL
                candidate_buffer[j] = candidate_point_stream.read();
                distance[j] = cal_dist_hw(q_temp_point, candidate_buffer[j]);
                index[j] = j;
            }
        }
        else
        {
            for(int j = 0; j < k_Octree_threshold; j++)
            {
#pragma HLS PIPELINE II=1
// #pragma HLS UNROLL
                distance[j] = cal_dist_hw(q_temp_point, candidate_buffer[j]); // 使用当前查询点计算距离
                index[j] = j;
            }
        }
        // bitonicSort128_hw(distance, index);
        // bitonic_sort_rec(distance, index, 0, k_max_point_num, true);
        // bitonicSort256_hw(distance, index);
        bitonicSort64_hw(distance, index);
        // quick_sort_hw(distance, index);
        // bitonicSort512_hw(distance, index);

        type_dist_hw candidate_res0 = distance[0];

        type_dist_hw candidate_res1 = distance[63];
        

        voxel_int nearest_index = (candidate_res0 < candidate_res1) ? index[0] : index[63];
        result[i] = candidate_buffer[nearest_index];

    }
}

void Ds_Oc_KNN_search_hw(My_PointXYZI_HW* result_hw, My_PointXYZI* query_set, voxel_int* q_set_start_index, int* q_flag,
    My_PointXYZI_HW* Oc_ordered_set, int* Oc_voxel_flag, int* Oc_voxel_first_index, int q_set_size,
    type_point_hw data_set_xmin_PL_hw, type_point_hw data_set_ymin_PL_hw, type_point_hw data_set_zmin_PL_hw, type_point_hw voxel_split_unit_PL_hw,
    voxel_int split_x_size_PL_hw, voxel_int split_y_size_PL_hw, voxel_int split_z_size_PL_hw, voxel_int total_voxel_size_PL_hw)
{
    #pragma HLS INTERFACE mode=m_axi        depth=k_query_set_size       port=result_hw              bundle=geme0
    #pragma HLS INTERFACE mode=m_axi        depth=k_query_set_size       port=query_set              bundle=geme1
    #pragma HLS INTERFACE mode=m_axi        depth=k_query_set_size       port=q_set_start_index      bundle=geme2
    #pragma HLS INTERFACE mode=m_axi        depth=k_query_set_size       port=q_flag                 bundle=geme3
    #pragma HLS INTERFACE mode=m_axi        depth=k_reference_set_size   port=Oc_ordered_set         bundle=geme4
    #pragma HLS INTERFACE mode=m_axi        depth=k_voxels_number_max*64 port=Oc_voxel_flag          bundle=geme5
    #pragma HLS INTERFACE mode=m_axi        depth=k_voxels_number_max*64 port=Oc_voxel_first_index   bundle=geme6

    #pragma HLS INTERFACE mode=s_axilite                    port=q_set_size             bundle=control
    #pragma HLS INTERFACE mode=s_axilite                    port=data_set_xmin_PL_hw    bundle=control
    #pragma HLS INTERFACE mode=s_axilite                    port=data_set_ymin_PL_hw    bundle=control
    #pragma HLS INTERFACE mode=s_axilite                    port=data_set_zmin_PL_hw    bundle=control
    #pragma HLS INTERFACE mode=s_axilite                    port=voxel_split_unit_PL_hw bundle=control
    #pragma HLS INTERFACE mode=s_axilite                    port=split_x_size_PL_hw     bundle=control
    #pragma HLS INTERFACE mode=s_axilite                    port=split_y_size_PL_hw     bundle=control
    #pragma HLS INTERFACE mode=s_axilite                    port=split_z_size_PL_hw     bundle=control
    #pragma HLS INTERFACE mode=s_axilite                    port=total_voxel_size_PL_hw bundle=control
    #pragma HLS INTERFACE mode=s_axilite                    port=return                 bundle=control



    split_x_array_size_hw = split_x_size_PL_hw;
    split_y_array_size_hw = split_y_size_PL_hw;
    split_z_array_size_hw = split_z_size_PL_hw;
    total_voxel_size_hw = total_voxel_size_PL_hw;
    data_set_xmin_hw = data_set_xmin_PL_hw;
    data_set_ymin_hw = data_set_ymin_PL_hw;
    data_set_zmin_hw = data_set_zmin_PL_hw;
    voxel_split_unit_hw = voxel_split_unit_PL_hw;

    hls::stream<My_PointXYZI_HW> q_set_stream;
    hls::stream<My_PointXYZI_HW> q_set_stream1;
    hls::stream<My_PointXYZI_HW> q_set_stream2;
    hls::stream<int> q_flag_stream;
    hls::stream<voxel_int> q_hash_stream;
    hls::stream<voxel_int> q_set_start_index_stream;
    hls::stream<voxel_int> q_set_start_index_stream1;
    hls::stream<My_PointXYZI_HW> candidate_point_stream;
    hls::stream<bool> update_flag_stream;

    #pragma HLS STREAM variable=q_set_stream depth=k_query_set_size
    #pragma HLS STREAM variable=q_set_stream1 depth=k_query_set_size
    #pragma HLS STREAM variable=q_set_stream2 depth=k_query_set_size
    #pragma HLS STREAM variable=q_flag_stream depth=k_query_set_size
    #pragma HLS STREAM variable=q_hash_stream depth=k_query_set_size
    #pragma HLS STREAM variable=q_set_start_index_stream depth=k_query_set_size
    #pragma HLS STREAM variable=q_set_start_index_stream1 depth=k_query_set_size
    #pragma HLS STREAM variable=candidate_point_stream depth=k_query_set_size * 20
    #pragma HLS STREAM variable=update_flag_stream depth=k_query_set_size

    #pragma HLS DATAFLOW

    // My_PointXYZI_HW candidate_point_hw[k_query_set_size * k_max_point_num];
    // #pragma HLS RESOURCE variable=candidate_point_hw core=RAM_2P_BRAM

    input_src_hw(query_set, q_set_stream, q_set_start_index, q_set_start_index_stream, q_flag, q_flag_stream, q_set_size);

    calculate_hash_hw(q_set_stream, q_set_stream1, q_hash_stream, q_set_size);

    cal_position_hw(q_set_stream1, q_set_stream2, q_hash_stream, q_flag_stream, q_set_start_index_stream, q_set_start_index_stream1, Oc_voxel_flag, Oc_voxel_first_index, q_set_size);

    cal_candidate_point_hw(candidate_point_stream, q_set_start_index_stream1, update_flag_stream, Oc_ordered_set, q_set_size);

    KNN_search_hw(result_hw, candidate_point_stream, update_flag_stream, q_set_stream2, q_set_size);

}
//quick_sort
void quick_sort_hw(type_dist_hw* arr, voxel_int* arr_index) {
    // 定义栈：存储待排序区间的左右边界
    const int MAX_STACK_SIZE = 64; // 最多存储64个区间
    int stack_left[MAX_STACK_SIZE];
    int stack_right[MAX_STACK_SIZE];

    int top = 0; // 栈顶指针

    // 初始化：将整个数组区间入栈
    stack_left[top] = 0;
    stack_right[top] = 63;
    top++;

    while (top > 0) {
        // 出栈
        top--;
        int low = stack_left[top];
        int high = stack_right[top];

        // 终止条件
        if (low >= high) continue;

        // 选择枢轴（中间元素）
        int pivot_idx = low + (high - low) / 2;
        type_dist_hw pivot_value = arr[pivot_idx];

        // 交换枢轴到末尾
        {
            type_dist_hw temp = arr[pivot_idx];
            arr[pivot_idx] = arr[high];
            arr[high] = temp;

            int temp_idx = arr_index[pivot_idx];
            arr_index[pivot_idx] = arr_index[high];
            arr_index[high] = temp_idx;
        }

        // 分区
        int i = low;
        for (int j = low; j < high; j++) {
#pragma HLS PIPELINE II=1
            if (arr[j] <= pivot_value) {
                // 交换 arr[i] 和 arr[j]
                {
                    type_dist_hw temp = arr[i];
                    arr[i] = arr[j];
                    arr[j] = temp;
                }
                // 交换索引
                {
                    int temp_idx = arr_index[i];
                    arr_index[i] = arr_index[j];
                    arr_index[j] = temp_idx;
                }
                i++;
            }
        }
        // 将枢轴放到正确位置
        {
            type_dist_hw temp = arr[i];
            arr[i] = arr[high];
            arr[high] = temp;
        }
        {
            int temp_idx = arr_index[i];
            arr_index[i] = arr_index[high];
            arr_index[high] = temp_idx;
        }

        // 入栈左区间
        if (low < i - 1) {
            stack_left[top] = low;
            stack_right[top] = i - 1;
            top++;
        }

        // 入栈右区间
        if (i + 1 < high) {
            stack_left[top] = i + 1;
            stack_right[top] = high;
            top++;
        }
    }
}
//bitonic_sort 512
void swap(type_dist_hw &a, type_dist_hw &b, voxel_int &a_idx, voxel_int &b_idx) {
#pragma HLS inline
    type_dist_hw temp = a;
    a = b;
    b = temp;

    voxel_int temp_idx = a_idx;
    a_idx = b_idx;
    b_idx = temp_idx;
}

void compare_swap(type_dist_hw &a, type_dist_hw &b, voxel_int &a_idx, voxel_int &b_idx, bool dir) {
#pragma HLS inline
    if (dir == (a > b)) {
        swap(a, b, a_idx, b_idx);
    }
}

void bitonic_merge(type_dist_hw arr[512], voxel_int arr_index[512], int low, int cnt, bool dir) {
#pragma HLS inline off
    if (cnt > 1) {
        int k = cnt / 2;
        for (int i = low; i < low + k; i++) {
#pragma HLS pipeline II=1
            compare_swap(arr[i], arr[i + k], arr_index[i], arr_index[i + k], dir);
        }
        bitonic_merge(arr, arr_index, low, k, dir);
        bitonic_merge(arr, arr_index, low + k, k, dir);
    }
}

void bitonic_sort_rec(type_dist_hw arr[512], voxel_int arr_index[512], int low, int cnt, bool dir) {
#pragma HLS inline off
    if (cnt > 1) {
        int k = cnt / 2;
        bitonic_sort_rec(arr, arr_index, low, k, true);      // 前半段升序
        bitonic_sort_rec(arr, arr_index, low + k, k, false); // 后半段降序
        bitonic_merge(arr, arr_index, low, cnt, dir);        // 合并整个序列成升序或降序
    }
}

//bitonic_sort
void compareSwap(type_dist_hw &a, type_dist_hw &b, voxel_int &a_index, voxel_int &b_index, bool dir) {
    #pragma HLS INLINE
    if (dir == (a > b)) {
        type_dist_hw tmp = a;
        a = b;
        b = tmp;

        voxel_int tmp_index = a_index;
        a_index = b_index;
        b_index = tmp_index;
    }
}

void bitonicSort512_hw(type_dist_hw *arr, voxel_int *arr_index) {
    // 迭代实现，不使用递归
    const int N = k_max_point_num;

    // Bitonic sort主循环，k为分块长度
    for (int k = 2; k <= N; k = k << 1) {
        // 子块比较距离
        for (int j = k >> 1; j > 0; j = j >> 1) {
            // 并行比较交换
            for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
                int ixj = i ^ j;
                if (ixj > i) {
                    // 升序段
                    if ((i & k) == 0) {
                        if (arr[i] > arr[ixj]) {
                            // 交换数据和索引
                            type_dist_hw tmp = arr[i];
                            arr[i] = arr[ixj];
                            arr[ixj] = tmp;
                            voxel_int tmp_idx = arr_index[i];
                            arr_index[i] = arr_index[ixj];
                            arr_index[ixj] = tmp_idx;
                        }
                    }
                    // 降序段
                    else {
                        if (arr[i] < arr[ixj]) {
                            // 交换数据和索引
                            type_dist_hw tmp = arr[i];
                            arr[i] = arr[ixj];
                            arr[ixj] = tmp;
                            voxel_int tmp_idx = arr_index[i];
                            arr_index[i] = arr_index[ixj];
                            arr_index[ixj] = tmp_idx;
                        }
                    }
                }
            }
        }
    }
}

//
void bitonicSort128_hw(type_dist_hw* arr, voxel_int* arr_index) 
{
// 对数组进行完全分区以实现并行访问
#pragma HLS ARRAY_PARTITION variable=arr complete
#pragma HLS ARRAY_PARTITION variable=arr_index complete
    //Stage 1
    for(int i = 0; i < 64; i++) {
#pragma HLS UNROLL factor = 16
        compareSwap(arr[i*2], arr[i*2 + 1], arr_index[i*2], arr_index[i*2 + 1], i%2==0);
    }

    //Stage 2
    for(int j = 0; j < 32; j++) {
#pragma HLS UNROLL factor = 16
        compareSwap(arr[j*4], arr[j*4 + 2], arr_index[j*4], arr_index[j*4 + 2], j%2==0);
        compareSwap(arr[j*4 + 1], arr[j*4 + 3], arr_index[j*4 + 1], arr_index[j*4 + 3], j%2==0);
    }
    for(int jj = 0; jj < 32; jj++) {
#pragma HLS UNROLL factor = 16
        compareSwap(arr[jj*4], arr[jj*4 + 1], arr_index[jj*4], arr_index[jj*4 + 1], jj%2==0);
        compareSwap(arr[jj*4 + 2], arr[jj*4 + 3], arr_index[jj*4 + 2], arr_index[jj*4 + 3], jj%2==0);
    }

    //Stage 3
    for(int k = 0; k < 16; k++) {
#pragma HLS UNROLL factor = 8
        compareSwap(arr[k*8], arr[k*8 + 4], arr_index[k*8], arr_index[k*8 + 4], k%2==0);
        compareSwap(arr[k*8 + 1], arr[k*8 + 5], arr_index[k*8 + 1], arr_index[k*8 + 5], k%2==0);
        compareSwap(arr[k*8 + 2], arr[k*8 + 6], arr_index[k*8 + 2], arr_index[k*8 + 6], k%2==0);
        compareSwap(arr[k*8 + 3], arr[k*8 + 7], arr_index[k*8 + 3], arr_index[k*8 + 7], k%2==0);
    }
    for(int kk = 0; kk < 16; kk++) {
#pragma HLS UNROLL factor = 8
        compareSwap(arr[kk*8], arr[kk*8 + 2], arr_index[kk*8], arr_index[kk*8 + 2], kk%2==0);
        compareSwap(arr[kk*8 + 1], arr[kk*8 + 3], arr_index[kk*8 + 1], arr_index[kk*8 + 3], kk%2==0);
        compareSwap(arr[kk*8 + 4], arr[kk*8 + 6], arr_index[kk*8 + 4], arr_index[kk*8 + 6], kk%2==0);
        compareSwap(arr[kk*8 + 5], arr[kk*8 + 7], arr_index[kk*8 + 5], arr_index[kk*8 + 7], kk%2==0);
    }
    for(int kkk = 0; kkk < 16; kkk++) {
#pragma HLS UNROLL factor = 8
        compareSwap(arr[kkk*8], arr[kkk*8 + 1], arr_index[kkk*8], arr_index[kkk*8 + 1], kkk%2==0);
        compareSwap(arr[kkk*8 + 2], arr[kkk*8 + 3], arr_index[kkk*8 + 2], arr_index[kkk*8 + 3], kkk%2==0);
        compareSwap(arr[kkk*8 + 4], arr[kkk*8 + 5], arr_index[kkk*8 + 4], arr_index[kkk*8 + 5], kkk%2==0);
        compareSwap(arr[kkk*8 + 6], arr[kkk*8 + 7], arr_index[kkk*8 + 6], arr_index[kkk*8 + 7], kkk%2==0);
    }

    //Stage 4
    for(int l = 0; l < 8; l++){
#pragma HLS UNROLL factor = 4
        compareSwap(arr[l*16], arr[l*16 + 8], arr_index[l*16], arr_index[l*16 + 8], l%2==0);
        compareSwap(arr[l*16 + 1], arr[l*16 + 9], arr_index[l*16 + 1], arr_index[l*16 + 9], l%2==0);
        compareSwap(arr[l*16 + 2], arr[l*16 + 10], arr_index[l*16 + 2], arr_index[l*16 + 10], l%2==0);
        compareSwap(arr[l*16 + 3], arr[l*16 + 11], arr_index[l*16 + 3], arr_index[l*16 + 11], l%2==0);
        compareSwap(arr[l*16 + 4], arr[l*16 + 12], arr_index[l*16 + 4], arr_index[l*16 + 12], l%2==0);
        compareSwap(arr[l*16 + 5], arr[l*16 + 13], arr_index[l*16 + 5], arr_index[l*16 + 13], l%2==0);
        compareSwap(arr[l*16 + 6], arr[l*16 + 14], arr_index[l*16 + 6], arr_index[l*16 + 14], l%2==0);
        compareSwap(arr[l*16 + 7], arr[l*16 + 15], arr_index[l*16 + 7], arr_index[l*16 + 15], l%2==0);
    }

    for(int ll = 0; ll < 8; ll++){
#pragma HLS UNROLL factor = 4
        compareSwap(arr[ll*16], arr[ll*16 + 4], arr_index[ll*16], arr_index[ll*16 + 4], ll%2==0);
        compareSwap(arr[ll*16 + 1], arr[ll*16 + 5], arr_index[ll*16 + 1], arr_index[ll*16 + 5], ll%2==0);
        compareSwap(arr[ll*16 + 2], arr[ll*16 + 6], arr_index[ll*16 + 2], arr_index[ll*16 + 6], ll%2==0);
        compareSwap(arr[ll*16 + 3], arr[ll*16 + 7], arr_index[ll*16 + 3], arr_index[ll*16 + 7], ll%2==0);
        compareSwap(arr[ll*16 + 8], arr[ll*16 + 12], arr_index[ll*16 + 8], arr_index[ll*16 + 12], ll%2==0);
        compareSwap(arr[ll*16 + 9], arr[ll*16 + 13], arr_index[ll*16 + 9], arr_index[ll*16 + 13], ll%2==0);
        compareSwap(arr[ll*16 + 10], arr[ll*16 + 14], arr_index[ll*16 + 10], arr_index[ll*16 + 14], ll%2==0);
        compareSwap(arr[ll*16 + 11], arr[ll*16 + 15], arr_index[ll*16 + 11], arr_index[ll*16 + 15], ll%2==0);
    }

    for (int lll = 0; lll < 8; lll++){
#pragma HLS UNROLL factor = 4
        compareSwap(arr[lll*16], arr[lll*16 + 2], arr_index[lll*16], arr_index[lll*16 + 2], lll%2==0);
        compareSwap(arr[lll*16 + 1], arr[lll*16 + 3], arr_index[lll*16 + 1], arr_index[lll*16 + 3], lll%2==0);
        compareSwap(arr[lll*16 + 4], arr[lll*16 + 6], arr_index[lll*16 + 4], arr_index[lll*16 + 6], lll%2==0);
        compareSwap(arr[lll*16 + 5], arr[lll*16 + 7], arr_index[lll*16 + 5], arr_index[lll*16 + 7], lll%2==0);
        compareSwap(arr[lll*16 + 8], arr[lll*16 + 10], arr_index[lll*16 + 8], arr_index[lll*16 + 10], lll%2==0);
        compareSwap(arr[lll*16 + 9], arr[lll*16 + 11], arr_index[lll*16 + 9], arr_index[lll*16 + 11], lll%2==0);
        compareSwap(arr[lll*16 + 12], arr[lll*16 + 14], arr_index[lll*16 + 12], arr_index[lll*16 + 14], lll%2==0);
        compareSwap(arr[lll*16 + 13], arr[lll*16 + 15], arr_index[lll*16 + 13], arr_index[lll*16 + 15], lll%2==0);
    }

    for (int llll = 0; llll < 8; llll++){
#pragma HLS UNROLL factor = 4
        compareSwap(arr[llll*16], arr[llll*16 + 1], arr_index[llll*16], arr_index[llll*16 + 1], llll%2==0);
        compareSwap(arr[llll*16 + 2], arr[llll*16 + 3], arr_index[llll*16 + 2], arr_index[llll*16 + 3], llll%2==0);
        compareSwap(arr[llll*16 + 4], arr[llll*16 + 5], arr_index[llll*16 + 4], arr_index[llll*16 + 5], llll%2==0);
        compareSwap(arr[llll*16 + 6], arr[llll*16 + 7], arr_index[llll*16 + 6], arr_index[llll*16 + 7], llll%2==0);
        compareSwap(arr[llll*16 + 8], arr[llll*16 + 9], arr_index[llll*16 + 8], arr_index[llll*16 + 9], llll%2==0);
        compareSwap(arr[llll*16 + 10], arr[llll*16 + 11], arr_index[llll*16 + 10], arr_index[llll*16 + 11], llll%2==0);
        compareSwap(arr[llll*16 + 12], arr[llll*16 + 13], arr_index[llll*16 + 12], arr_index[llll*16 + 13], llll%2==0);
        compareSwap(arr[llll*16 + 14], arr[llll*16 + 15], arr_index[llll*16 + 14], arr_index[llll*16 + 15], llll%2==0);
    }

    //Stage 5

    for(int m = 0; m < 4; m++){
#pragma HLS UNROLL
        compareSwap(arr[m*32], arr[m*32 + 16], arr_index[m*32], arr_index[m*32 + 16], m%2==0);
        compareSwap(arr[m*32 + 1], arr[m*32 + 17], arr_index[m*32 + 1], arr_index[m*32 + 17], m%2==0);
        compareSwap(arr[m*32 + 2], arr[m*32 + 18], arr_index[m*32 + 2], arr_index[m*32 + 18], m%2==0);
        compareSwap(arr[m*32 + 3], arr[m*32 + 19], arr_index[m*32 + 3], arr_index[m*32 + 19], m%2==0);
        compareSwap(arr[m*32 + 4], arr[m*32 + 20], arr_index[m*32 + 4], arr_index[m*32 + 20], m%2==0);
        compareSwap(arr[m*32 + 5], arr[m*32 + 21], arr_index[m*32 + 5], arr_index[m*32 + 21], m%2==0);
        compareSwap(arr[m*32 + 6], arr[m*32 + 22], arr_index[m*32 + 6], arr_index[m*32 + 22], m%2==0);
        compareSwap(arr[m*32 + 7], arr[m*32 + 23], arr_index[m*32 + 7], arr_index[m*32 + 23], m%2==0);
        compareSwap(arr[m*32 + 8], arr[m*32 + 24], arr_index[m*32 + 8], arr_index[m*32 + 24], m%2==0);
        compareSwap(arr[m*32 + 9], arr[m*32 + 25], arr_index[m*32 + 9], arr_index[m*32 + 25], m%2==0);
        compareSwap(arr[m*32 + 10], arr[m*32 + 26], arr_index[m*32 + 10], arr_index[m*32 + 26], m%2==0);
        compareSwap(arr[m*32 + 11], arr[m*32 + 27], arr_index[m*32 + 11], arr_index[m*32 + 27], m%2==0);
        compareSwap(arr[m*32 + 12], arr[m*32 + 28], arr_index[m*32 + 12], arr_index[m*32 + 28], m%2==0);
        compareSwap(arr[m*32 + 13], arr[m*32 + 29], arr_index[m*32 + 13], arr_index[m*32 + 29], m%2==0);
        compareSwap(arr[m*32 + 14], arr[m*32 + 30], arr_index[m*32 + 14], arr_index[m*32 + 30], m%2==0);
        compareSwap(arr[m*32 + 15], arr[m*32 + 31], arr_index[m*32 + 15], arr_index[m*32 + 31], m%2==0);
    }

    for(int mm = 0; mm < 4; mm++){
#pragma HLS UNROLL
        compareSwap(arr[mm*32], arr[mm*32 + 8], arr_index[mm*32], arr_index[mm*32 + 8], mm%2==0);
        compareSwap(arr[mm*32 + 1], arr[mm*32 + 9], arr_index[mm*32 + 1], arr_index[mm*32 + 9], mm%2==0);
        compareSwap(arr[mm*32 + 2], arr[mm*32 + 10], arr_index[mm*32 + 2], arr_index[mm*32 + 10], mm%2==0);
        compareSwap(arr[mm*32 + 3], arr[mm*32 + 11], arr_index[mm*32 + 3], arr_index[mm*32 + 11], mm%2==0);
        compareSwap(arr[mm*32 + 4], arr[mm*32 + 12], arr_index[mm*32 + 4], arr_index[mm*32 + 12], mm%2==0);
        compareSwap(arr[mm*32 + 5], arr[mm*32 + 13], arr_index[mm*32 + 5], arr_index[mm*32 + 13], mm%2==0);
        compareSwap(arr[mm*32 + 6], arr[mm*32 + 14], arr_index[mm*32 + 6], arr_index[mm*32 + 14], mm%2==0);
        compareSwap(arr[mm*32 + 7], arr[mm*32 + 15], arr_index[mm*32 + 7], arr_index[mm*32 + 15], mm%2==0);
        compareSwap(arr[mm*32 + 16], arr[mm*32 + 24], arr_index[mm*32 + 16], arr_index[mm*32 + 24], mm%2==0);
        compareSwap(arr[mm*32 + 17], arr[mm*32 + 25], arr_index[mm*32 + 17], arr_index[mm*32 + 25], mm%2==0);
        compareSwap(arr[mm*32 + 18], arr[mm*32 + 26], arr_index[mm*32 + 18], arr_index[mm*32 + 26], mm%2==0);
        compareSwap(arr[mm*32 + 19], arr[mm*32 + 27], arr_index[mm*32 + 19], arr_index[mm*32 + 27], mm%2==0);
        compareSwap(arr[mm*32 + 20], arr[mm*32 + 28], arr_index[mm*32 + 20], arr_index[mm*32 + 28], mm%2==0);
        compareSwap(arr[mm*32 + 21], arr[mm*32 + 29], arr_index[mm*32 + 21], arr_index[mm*32 + 29], mm%2==0);
        compareSwap(arr[mm*32 + 22], arr[mm*32 + 30], arr_index[mm*32 + 22], arr_index[mm*32 + 30], mm%2==0);
        compareSwap(arr[mm*32 + 23], arr[mm*32 + 31], arr_index[mm*32 + 23], arr_index[mm*32 + 31], mm%2==0);
    }

    for(int mmm = 0; mmm < 4; mmm++){
#pragma HLS UNROLL
        compareSwap(arr[mmm*32], arr[mmm*32 + 4], arr_index[mmm*32], arr_index[mmm*32 + 4], mmm%2==0);
        compareSwap(arr[mmm*32 + 1], arr[mmm*32 + 5], arr_index[mmm*32 + 1], arr_index[mmm*32 + 5], mmm%2==0);
        compareSwap(arr[mmm*32 + 2], arr[mmm*32 + 6], arr_index[mmm*32 + 2], arr_index[mmm*32 + 6], mmm%2==0);
        compareSwap(arr[mmm*32 + 3], arr[mmm*32 + 7], arr_index[mmm*32 + 3], arr_index[mmm*32 + 7], mmm%2==0);
        compareSwap(arr[mmm*32 + 8], arr[mmm*32 + 12], arr_index[mmm*32 + 8], arr_index[mmm*32 + 12], mmm%2==0);
        compareSwap(arr[mmm*32 + 9], arr[mmm*32 + 13], arr_index[mmm*32 + 9], arr_index[mmm*32 + 13], mmm%2==0);
        compareSwap(arr[mmm*32 + 10], arr[mmm*32 + 14], arr_index[mmm*32 + 10], arr_index[mmm*32 + 14], mmm%2==0);
        compareSwap(arr[mmm*32 + 11], arr[mmm*32 + 15], arr_index[mmm*32 + 11], arr_index[mmm*32 + 15], mmm%2==0);
        compareSwap(arr[mmm*32 + 16], arr[mmm*32 + 20], arr_index[mmm*32 + 16], arr_index[mmm*32 + 20], mmm%2==0);
        compareSwap(arr[mmm*32 + 17], arr[mmm*32 + 21], arr_index[mmm*32 + 17], arr_index[mmm*32 + 21], mmm%2==0);
        compareSwap(arr[mmm*32 + 18], arr[mmm*32 + 22], arr_index[mmm*32 + 18], arr_index[mmm*32 + 22], mmm%2==0);
        compareSwap(arr[mmm*32 + 19], arr[mmm*32 + 23], arr_index[mmm*32 + 19], arr_index[mmm*32 + 23], mmm%2==0);
        compareSwap(arr[mmm*32 + 24], arr[mmm*32 + 28], arr_index[mmm*32 + 24], arr_index[mmm*32 + 28], mmm%2==0);
        compareSwap(arr[mmm*32 + 25], arr[mmm*32 + 29], arr_index[mmm*32 + 25], arr_index[mmm*32 + 29], mmm%2==0);
        compareSwap(arr[mmm*32 + 26], arr[mmm*32 + 30], arr_index[mmm*32 + 26], arr_index[mmm*32 + 30], mmm%2==0);
        compareSwap(arr[mmm*32 + 27], arr[mmm*32 + 31], arr_index[mmm*32 + 27], arr_index[mmm*32 + 31], mmm%2==0);
    }

    for(int mmmm = 0; mmmm < 4; mmmm ++){
#pragma HLS UNROLL
        compareSwap(arr[mmmm*32], arr[mmmm*32 + 2], arr_index[mmmm*32], arr_index[mmmm*32 + 2], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 1], arr[mmmm*32 + 3], arr_index[mmmm*32 + 1], arr_index[mmmm*32 + 3], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 4], arr[mmmm*32 + 6], arr_index[mmmm*32 + 4], arr_index[mmmm*32 + 6], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 5], arr[mmmm*32 + 7], arr_index[mmmm*32 + 5], arr_index[mmmm*32 + 7], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 8], arr[mmmm*32 + 10], arr_index[mmmm*32 + 8], arr_index[mmmm*32 + 10], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 9], arr[mmmm*32 + 11], arr_index[mmmm*32 + 9], arr_index[mmmm*32 + 11], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 12], arr[mmmm*32 + 14], arr_index[mmmm*32 + 12], arr_index[mmmm*32 + 14], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 13], arr[mmmm*32 + 15], arr_index[mmmm*32 + 13], arr_index[mmmm*32 + 15], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 16], arr[mmmm*32 + 18], arr_index[mmmm*32 + 16], arr_index[mmmm*32 + 18], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 17], arr[mmmm*32 + 19], arr_index[mmmm*32 + 17], arr_index[mmmm*32 + 19], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 20], arr[mmmm*32 + 22], arr_index[mmmm*32 + 20], arr_index[mmmm*32 + 22], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 21], arr[mmmm*32 + 23], arr_index[mmmm*32 + 21], arr_index[mmmm*32 + 23], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 24], arr[mmmm*32 + 26], arr_index[mmmm*32 + 24], arr_index[mmmm*32 + 26], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 25], arr[mmmm*32 + 27], arr_index[mmmm*32 + 25], arr_index[mmmm*32 + 27], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 28], arr[mmmm*32 + 30], arr_index[mmmm*32 + 28], arr_index[mmmm*32 + 30], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 29], arr[mmmm*32 + 31], arr_index[mmmm*32 + 29], arr_index[mmmm*32 + 31], mmmm%2==0);
    }

    for(int m5 = 0; m5 < 4; m5++){
#pragma HLS UNROLL
        compareSwap(arr[m5*32], arr[m5*32 + 1], arr_index[m5*32], arr_index[m5*32 + 1], m5%2==0);
        compareSwap(arr[m5*32 + 2], arr[m5*32 + 3], arr_index[m5*32 + 2], arr_index[m5*32 + 3], m5%2==0);
        compareSwap(arr[m5*32 + 4], arr[m5*32 + 5], arr_index[m5*32 + 4], arr_index[m5*32 + 5], m5%2==0);
        compareSwap(arr[m5*32 + 6], arr[m5*32 + 7], arr_index[m5*32 + 6], arr_index[m5*32 + 7], m5%2==0);
        compareSwap(arr[m5*32 + 8], arr[m5*32 + 9], arr_index[m5*32 + 8], arr_index[m5*32 + 9], m5%2==0);
        compareSwap(arr[m5*32 + 10], arr[m5*32 + 11], arr_index[m5*32 + 10], arr_index[m5*32 + 11], m5%2==0);
        compareSwap(arr[m5*32 + 12], arr[m5*32 + 13], arr_index[m5*32 + 12], arr_index[m5*32 + 13], m5%2==0);
        compareSwap(arr[m5*32 + 14], arr[m5*32 + 15], arr_index[m5*32 + 14], arr_index[m5*32 + 15], m5%2==0);
        compareSwap(arr[m5*32 + 16], arr[m5*32 + 17], arr_index[m5*32 + 16], arr_index[m5*32 + 17], m5%2==0);
        compareSwap(arr[m5*32 + 18], arr[m5*32 + 19], arr_index[m5*32 + 18], arr_index[m5*32 + 19], m5%2==0);
        compareSwap(arr[m5*32 + 20], arr[m5*32 + 21], arr_index[m5*32 + 20], arr_index[m5*32 + 21], m5%2==0);
        compareSwap(arr[m5*32 + 22], arr[m5*32 + 23], arr_index[m5*32 + 22], arr_index[m5*32 + 23], m5%2==0);
        compareSwap(arr[m5*32 + 24], arr[m5*32 + 25], arr_index[m5*32 + 24], arr_index[m5*32 + 25], m5%2==0);
        compareSwap(arr[m5*32 + 26], arr[m5*32 + 27], arr_index[m5*32 + 26], arr_index[m5*32 + 27], m5%2==0);
        compareSwap(arr[m5*32 + 28], arr[m5*32 + 29], arr_index[m5*32 + 28], arr_index[m5*32 + 29], m5%2==0);
        compareSwap(arr[m5*32 + 30], arr[m5*32 + 31], arr_index[m5*32 + 30], arr_index[m5*32 + 31], m5%2==0);
    }

    //Stage 6

    for(int n = 0; n < 2; n++){

        compareSwap(arr[n*64], arr[n*64 + 32], arr_index[n*64], arr_index[n*64 + 32], n%2==0);
        compareSwap(arr[n*64 + 1], arr[n*64 + 33], arr_index[n*64 + 1], arr_index[n*64 + 33], n%2==0);
        compareSwap(arr[n*64 + 2], arr[n*64 + 34], arr_index[n*64 + 2], arr_index[n*64 + 34], n%2==0);
        compareSwap(arr[n*64 + 3], arr[n*64 + 35], arr_index[n*64 + 3], arr_index[n*64 + 35], n%2==0);
        compareSwap(arr[n*64 + 4], arr[n*64 + 36], arr_index[n*64 + 4], arr_index[n*64 + 36], n%2==0);
        compareSwap(arr[n*64 + 5], arr[n*64 + 37], arr_index[n*64 + 5], arr_index[n*64 + 37], n%2==0);
        compareSwap(arr[n*64 + 6], arr[n*64 + 38], arr_index[n*64 + 6], arr_index[n*64 + 38], n%2==0);
        compareSwap(arr[n*64 + 7], arr[n*64 + 39], arr_index[n*64 + 7], arr_index[n*64 + 39], n%2==0);
        compareSwap(arr[n*64 + 8], arr[n*64 + 40], arr_index[n*64 + 8], arr_index[n*64 + 40], n%2==0);
        compareSwap(arr[n*64 + 9], arr[n*64 + 41], arr_index[n*64 + 9], arr_index[n*64 + 41], n%2==0);
        compareSwap(arr[n*64 + 10], arr[n*64 + 42], arr_index[n*64 + 10], arr_index[n*64 + 42], n%2==0);
        compareSwap(arr[n*64 + 11], arr[n*64 + 43], arr_index[n*64 + 11], arr_index[n*64 + 43], n%2==0);
        compareSwap(arr[n*64 + 12], arr[n*64 + 44], arr_index[n*64 + 12], arr_index[n*64 + 44], n%2==0);
        compareSwap(arr[n*64 + 13], arr[n*64 + 45], arr_index[n*64 + 13], arr_index[n*64 + 45], n%2==0);
        compareSwap(arr[n*64 + 14], arr[n*64 + 46], arr_index[n*64 + 14], arr_index[n*64 + 46], n%2==0);
        compareSwap(arr[n*64 + 15], arr[n*64 + 47], arr_index[n*64 + 15], arr_index[n*64 + 47], n%2==0);
        compareSwap(arr[n*64 + 16], arr[n*64 + 48], arr_index[n*64 + 16], arr_index[n*64 + 48], n%2==0);
        compareSwap(arr[n*64 + 17], arr[n*64 + 49], arr_index[n*64 + 17], arr_index[n*64 + 49], n%2==0);
        compareSwap(arr[n*64 + 18], arr[n*64 + 50], arr_index[n*64 + 18], arr_index[n*64 + 50], n%2==0);
        compareSwap(arr[n*64 + 19], arr[n*64 + 51], arr_index[n*64 + 19], arr_index[n*64 + 51], n%2==0);
        compareSwap(arr[n*64 + 20], arr[n*64 + 52], arr_index[n*64 + 20], arr_index[n*64 + 52], n%2==0);
        compareSwap(arr[n*64 + 21], arr[n*64 + 53], arr_index[n*64 + 21], arr_index[n*64 + 53], n%2==0);
        compareSwap(arr[n*64 + 22], arr[n*64 + 54], arr_index[n*64 + 22], arr_index[n*64 + 54], n%2==0);
        compareSwap(arr[n*64 + 23], arr[n*64 + 55], arr_index[n*64 + 23], arr_index[n*64 + 55], n%2==0);
        compareSwap(arr[n*64 + 24], arr[n*64 + 56], arr_index[n*64 + 24], arr_index[n*64 + 56], n%2==0);
        compareSwap(arr[n*64 + 25], arr[n*64 + 57], arr_index[n*64 + 25], arr_index[n*64 + 57], n%2==0);
        compareSwap(arr[n*64 + 26], arr[n*64 + 58], arr_index[n*64 + 26], arr_index[n*64 + 58], n%2==0);
        compareSwap(arr[n*64 + 27], arr[n*64 + 59], arr_index[n*64 + 27], arr_index[n*64 + 59], n%2==0);
        compareSwap(arr[n*64 + 28], arr[n*64 + 60], arr_index[n*64 + 28], arr_index[n*64 + 60], n%2==0);
        compareSwap(arr[n*64 + 29], arr[n*64 + 61], arr_index[n*64 + 29], arr_index[n*64 + 61], n%2==0);
        compareSwap(arr[n*64 + 30], arr[n*64 + 62], arr_index[n*64 + 30], arr_index[n*64 + 62], n%2==0);
        compareSwap(arr[n*64 + 31], arr[n*64 + 63], arr_index[n*64 + 31], arr_index[n*64 + 63], n%2==0);
    }

    for(int n2 = 0; n2 < 2; n2++){

        compareSwap(arr[n2*64], arr[n2*64 + 16], arr_index[n2*64], arr_index[n2*64 + 16], n2%2==0);
        compareSwap(arr[n2*64 + 1], arr[n2*64 + 17], arr_index[n2*64 + 1], arr_index[n2*64 + 17], n2%2==0);
        compareSwap(arr[n2*64 + 2], arr[n2*64 + 18], arr_index[n2*64 + 2], arr_index[n2*64 + 18], n2%2==0);
        compareSwap(arr[n2*64 + 3], arr[n2*64 + 19], arr_index[n2*64 + 3], arr_index[n2*64 + 19], n2%2==0);
        compareSwap(arr[n2*64 + 4], arr[n2*64 + 20], arr_index[n2*64 + 4], arr_index[n2*64 + 20], n2%2==0);
        compareSwap(arr[n2*64 + 5], arr[n2*64 + 21], arr_index[n2*64 + 5], arr_index[n2*64 + 21], n2%2==0);
        compareSwap(arr[n2*64 + 6], arr[n2*64 + 22], arr_index[n2*64 + 6], arr_index[n2*64 + 22], n2%2==0);
        compareSwap(arr[n2*64 + 7], arr[n2*64 + 23], arr_index[n2*64 + 7], arr_index[n2*64 + 23], n2%2==0);
        compareSwap(arr[n2*64 + 8], arr[n2*64 + 24], arr_index[n2*64 + 8], arr_index[n2*64 + 24], n2%2==0);
        compareSwap(arr[n2*64 + 9], arr[n2*64 + 25], arr_index[n2*64 + 9], arr_index[n2*64 + 25], n2%2==0);
        compareSwap(arr[n2*64 + 10], arr[n2*64 + 26], arr_index[n2*64 + 10], arr_index[n2*64 + 26], n2%2==0);
        compareSwap(arr[n2*64 + 11], arr[n2*64 + 27], arr_index[n2*64 + 11], arr_index[n2*64 + 27], n2%2==0);
        compareSwap(arr[n2*64 + 12], arr[n2*64 + 28], arr_index[n2*64 + 12], arr_index[n2*64 + 28], n2%2==0);
        compareSwap(arr[n2*64 + 13], arr[n2*64 + 29], arr_index[n2*64 + 13], arr_index[n2*64 + 29], n2%2==0);
        compareSwap(arr[n2*64 + 14], arr[n2*64 + 30], arr_index[n2*64 + 14], arr_index[n2*64 + 30], n2%2==0);
        compareSwap(arr[n2*64 + 15], arr[n2*64 + 31], arr_index[n2*64 + 15], arr_index[n2*64 + 31], n2%2==0);
        
        compareSwap(arr[n2*64 + 32], arr[n2*64 + 48], arr_index[n2*64 + 32], arr_index[n2*64 + 48], n2%2==0);
        compareSwap(arr[n2*64 + 33], arr[n2*64 + 49], arr_index[n2*64 + 33], arr_index[n2*64 + 49], n2%2==0);
        compareSwap(arr[n2*64 + 34], arr[n2*64 + 50], arr_index[n2*64 + 34], arr_index[n2*64 + 50], n2%2==0);
        compareSwap(arr[n2*64 + 35], arr[n2*64 + 51], arr_index[n2*64 + 35], arr_index[n2*64 + 51], n2%2==0);
        compareSwap(arr[n2*64 + 36], arr[n2*64 + 52], arr_index[n2*64 + 36], arr_index[n2*64 + 52], n2%2==0);
        compareSwap(arr[n2*64 + 37], arr[n2*64 + 53], arr_index[n2*64 + 37], arr_index[n2*64 + 53], n2%2==0);
        compareSwap(arr[n2*64 + 38], arr[n2*64 + 54], arr_index[n2*64 + 38], arr_index[n2*64 + 54], n2%2==0);
        compareSwap(arr[n2*64 + 39], arr[n2*64 + 55], arr_index[n2*64 + 39], arr_index[n2*64 + 55], n2%2==0);
        compareSwap(arr[n2*64 + 40], arr[n2*64 + 56], arr_index[n2*64 + 40], arr_index[n2*64 + 56], n2%2==0);
        compareSwap(arr[n2*64 + 41], arr[n2*64 + 57], arr_index[n2*64 + 41], arr_index[n2*64 + 57], n2%2==0);
        compareSwap(arr[n2*64 + 42], arr[n2*64 + 58], arr_index[n2*64 + 42], arr_index[n2*64 + 58], n2%2==0);
        compareSwap(arr[n2*64 + 43], arr[n2*64 + 59], arr_index[n2*64 + 43], arr_index[n2*64 + 59], n2%2==0);
        compareSwap(arr[n2*64 + 44], arr[n2*64 + 60], arr_index[n2*64 + 44], arr_index[n2*64 + 60], n2%2==0);
        compareSwap(arr[n2*64 + 45], arr[n2*64 + 61], arr_index[n2*64 + 45], arr_index[n2*64 + 61], n2%2==0);
        compareSwap(arr[n2*64 + 46], arr[n2*64 + 62], arr_index[n2*64 + 46], arr_index[n2*64 + 62], n2%2==0);
        compareSwap(arr[n2*64 + 47], arr[n2*64 + 63], arr_index[n2*64 + 47], arr_index[n2*64 + 63], n2%2==0);
    }

    for(int n3 = 0; n3 < 2; n3++){

        compareSwap(arr[n3*64], arr[n3*64 + 8], arr_index[n3*64], arr_index[n3*64 + 8], n3%2==0);
        compareSwap(arr[n3*64 + 1], arr[n3*64 + 9], arr_index[n3*64 + 1], arr_index[n3*64 + 9], n3%2==0);
        compareSwap(arr[n3*64 + 2], arr[n3*64 + 10], arr_index[n3*64 + 2], arr_index[n3*64 + 10], n3%2==0);
        compareSwap(arr[n3*64 + 3], arr[n3*64 + 11], arr_index[n3*64 + 3], arr_index[n3*64 + 11], n3%2==0);
        compareSwap(arr[n3*64 + 4], arr[n3*64 + 12], arr_index[n3*64 + 4], arr_index[n3*64 + 12], n3%2==0);
        compareSwap(arr[n3*64 + 5], arr[n3*64 + 13], arr_index[n3*64 + 5], arr_index[n3*64 + 13], n3%2==0);
        compareSwap(arr[n3*64 + 6], arr[n3*64 + 14], arr_index[n3*64 + 6], arr_index[n3*64 + 14], n3%2==0);
        compareSwap(arr[n3*64 + 7], arr[n3*64 + 15], arr_index[n3*64 + 7], arr_index[n3*64 + 15], n3%2==0);

        compareSwap(arr[n3*64 + 16], arr[n3*64 + 24], arr_index[n3*64 + 16], arr_index[n3*64 + 24], n3%2==0);
        compareSwap(arr[n3*64 + 17], arr[n3*64 + 25], arr_index[n3*64 + 17], arr_index[n3*64 + 25], n3%2==0);
        compareSwap(arr[n3*64 + 18], arr[n3*64 + 26], arr_index[n3*64 + 18], arr_index[n3*64 + 26], n3%2==0);
        compareSwap(arr[n3*64 + 19], arr[n3*64 + 27], arr_index[n3*64 + 19], arr_index[n3*64 + 27], n3%2==0);
        compareSwap(arr[n3*64 + 20], arr[n3*64 + 28], arr_index[n3*64 + 20], arr_index[n3*64 + 28], n3%2==0);
        compareSwap(arr[n3*64 + 21], arr[n3*64 + 29], arr_index[n3*64 + 21], arr_index[n3*64 + 29], n3%2==0);
        compareSwap(arr[n3*64 + 22], arr[n3*64 + 30], arr_index[n3*64 + 22], arr_index[n3*64 + 30], n3%2==0);
        compareSwap(arr[n3*64 + 23], arr[n3*64 + 31], arr_index[n3*64 + 23], arr_index[n3*64 + 31], n3%2==0);

        compareSwap(arr[n3*64 + 32], arr[n3*64 + 40], arr_index[n3*64 + 32], arr_index[n3*64 + 40], n3%2==0);
        compareSwap(arr[n3*64 + 33], arr[n3*64 + 41], arr_index[n3*64 + 33], arr_index[n3*64 + 41], n3%2==0);
        compareSwap(arr[n3*64 + 34], arr[n3*64 + 42], arr_index[n3*64 + 34], arr_index[n3*64 + 42], n3%2==0);
        compareSwap(arr[n3*64 + 35], arr[n3*64 + 43], arr_index[n3*64 + 35], arr_index[n3*64 + 43], n3%2==0);
        compareSwap(arr[n3*64 + 36], arr[n3*64 + 44], arr_index[n3*64 + 36], arr_index[n3*64 + 44], n3%2==0);
        compareSwap(arr[n3*64 + 37], arr[n3*64 + 45], arr_index[n3*64 + 37], arr_index[n3*64 + 45], n3%2==0);
        compareSwap(arr[n3*64 + 38], arr[n3*64 + 46], arr_index[n3*64 + 38], arr_index[n3*64 + 46], n3%2==0);
        compareSwap(arr[n3*64 + 39], arr[n3*64 + 47], arr_index[n3*64 + 39], arr_index[n3*64 + 47], n3%2==0);
        
        compareSwap(arr[n3*64 + 48], arr[n3*64 + 56], arr_index[n3*64 + 48], arr_index[n3*64 + 56], n3%2==0);
        compareSwap(arr[n3*64 + 49], arr[n3*64 + 57], arr_index[n3*64 + 49], arr_index[n3*64 + 57], n3%2==0);
        compareSwap(arr[n3*64 + 50], arr[n3*64 + 58], arr_index[n3*64 + 50], arr_index[n3*64 + 58], n3%2==0);
        compareSwap(arr[n3*64 + 51], arr[n3*64 + 59], arr_index[n3*64 + 51], arr_index[n3*64 + 59], n3%2==0);
        compareSwap(arr[n3*64 + 52], arr[n3*64 + 60], arr_index[n3*64 + 52], arr_index[n3*64 + 60], n3%2==0);
        compareSwap(arr[n3*64 + 53], arr[n3*64 + 61], arr_index[n3*64 + 53], arr_index[n3*64 + 61], n3%2==0);
        compareSwap(arr[n3*64 + 54], arr[n3*64 + 62], arr_index[n3*64 + 54], arr_index[n3*64 + 62], n3%2==0);
        compareSwap(arr[n3*64 + 55], arr[n3*64 + 63], arr_index[n3*64 + 55], arr_index[n3*64 + 63], n3%2==0);
    }

    for(int n4 = 0; n4 < 2; n4++){

        compareSwap(arr[n4*64], arr[n4*64 + 4], arr_index[n4*64], arr_index[n4*64 + 4], n4%2==0);
        compareSwap(arr[n4*64 + 1], arr[n4*64 + 5], arr_index[n4*64 + 1], arr_index[n4*64 + 5], n4%2==0);
        compareSwap(arr[n4*64 + 2], arr[n4*64 + 6], arr_index[n4*64 + 2], arr_index[n4*64 + 6], n4%2==0);
        compareSwap(arr[n4*64 + 3], arr[n4*64 + 7], arr_index[n4*64 + 3], arr_index[n4*64 + 7], n4%2==0);
        compareSwap(arr[n4*64 + 8], arr[n4*64 + 12], arr_index[n4*64 + 8], arr_index[n4*64 + 12], n4%2==0);
        compareSwap(arr[n4*64 + 9], arr[n4*64 + 13], arr_index[n4*64 + 9], arr_index[n4*64 + 13], n4%2==0);
        compareSwap(arr[n4*64 + 10], arr[n4*64 + 14], arr_index[n4*64 + 10], arr_index[n4*64 + 14], n4%2==0);
        compareSwap(arr[n4*64 + 11], arr[n4*64 + 15], arr_index[n4*64 + 11], arr_index[n4*64 + 15], n4%2==0);
        compareSwap(arr[n4*64 + 16], arr[n4*64 + 20], arr_index[n4*64 + 16], arr_index[n4*64 + 20], n4%2==0);
        compareSwap(arr[n4*64 + 17], arr[n4*64 + 21], arr_index[n4*64 + 17], arr_index[n4*64 + 21], n4%2==0);
        compareSwap(arr[n4*64 + 18], arr[n4*64 + 22], arr_index[n4*64 + 18], arr_index[n4*64 + 22], n4%2==0);
        compareSwap(arr[n4*64 + 19], arr[n4*64 + 23], arr_index[n4*64 + 19], arr_index[n4*64 + 23], n4%2==0);
        compareSwap(arr[n4*64 + 24], arr[n4*64 + 28], arr_index[n4*64 + 24], arr_index[n4*64 + 28], n4%2==0);
        compareSwap(arr[n4*64 + 25], arr[n4*64 + 29], arr_index[n4*64 + 25], arr_index[n4*64 + 29], n4%2==0);
        compareSwap(arr[n4*64 + 26], arr[n4*64 + 30], arr_index[n4*64 + 26], arr_index[n4*64 + 30], n4%2==0);
        compareSwap(arr[n4*64 + 27], arr[n4*64 + 31], arr_index[n4*64 + 27], arr_index[n4*64 + 31], n4%2==0);
        compareSwap(arr[n4*64 + 32], arr[n4*64 + 36], arr_index[n4*64 + 32], arr_index[n4*64 + 36], n4%2==0);
        compareSwap(arr[n4*64 + 33], arr[n4*64 + 37], arr_index[n4*64 + 33], arr_index[n4*64 + 37], n4%2==0);
        compareSwap(arr[n4*64 + 34], arr[n4*64 + 38], arr_index[n4*64 + 34], arr_index[n4*64 + 38], n4%2==0);
        compareSwap(arr[n4*64 + 35], arr[n4*64 + 39], arr_index[n4*64 + 35], arr_index[n4*64 + 39], n4%2==0);

        compareSwap(arr[n4*64 + 40], arr[n4*64 + 44], arr_index[n4*64 + 40], arr_index[n4*64 + 44], n4%2==0);
        compareSwap(arr[n4*64 + 41], arr[n4*64 + 45], arr_index[n4*64 + 41], arr_index[n4*64 + 45], n4%2==0);
        compareSwap(arr[n4*64 + 42], arr[n4*64 + 46], arr_index[n4*64 + 42], arr_index[n4*64 + 46], n4%2==0);
        compareSwap(arr[n4*64 + 43], arr[n4*64 + 47], arr_index[n4*64 + 43], arr_index[n4*64 + 47], n4%2==0);

        compareSwap(arr[n4*64 + 48], arr[n4*64 + 52], arr_index[n4*64 + 48], arr_index[n4*64 + 52], n4%2==0);
        compareSwap(arr[n4*64 + 49], arr[n4*64 + 53], arr_index[n4*64 + 49], arr_index[n4*64 + 53], n4%2==0);
        compareSwap(arr[n4*64 + 50], arr[n4*64 + 54], arr_index[n4*64 + 50], arr_index[n4*64 + 54], n4%2==0);
        compareSwap(arr[n4*64 + 51], arr[n4*64 + 55], arr_index[n4*64 + 51], arr_index[n4*64 + 55], n4%2==0);

        compareSwap(arr[n4*64 + 56], arr[n4*64 + 60], arr_index[n4*64 + 56], arr_index[n4*64 + 60], n4%2==0);
        compareSwap(arr[n4*64 + 57], arr[n4*64 + 61], arr_index[n4*64 + 57], arr_index[n4*64 + 61], n4%2==0);
        compareSwap(arr[n4*64 + 58], arr[n4*64 + 62], arr_index[n4*64 + 58], arr_index[n4*64 + 62], n4%2==0);
        compareSwap(arr[n4*64 + 59], arr[n4*64 + 63], arr_index[n4*64 + 59], arr_index[n4*64 + 63], n4%2==0);
    }

    for(int n5 = 0; n5 < 2; n5++){

        compareSwap(arr[n5*64], arr[n5*64 + 2], arr_index[n5*64], arr_index[n5*64 + 2], n5%2==0);
        compareSwap(arr[n5*64 + 1], arr[n5*64 + 3], arr_index[n5*64 + 1], arr_index[n5*64 + 3], n5%2==0);
        compareSwap(arr[n5*64 + 4], arr[n5*64 + 6], arr_index[n5*64 + 4], arr_index[n5*64 + 6], n5%2==0);
        compareSwap(arr[n5*64 + 5], arr[n5*64 + 7], arr_index[n5*64 + 5], arr_index[n5*64 + 7], n5%2==0);
        compareSwap(arr[n5*64 + 8], arr[n5*64 + 10], arr_index[n5*64 + 8], arr_index[n5*64 + 10], n5%2==0);
        compareSwap(arr[n5*64 + 9], arr[n5*64 + 11], arr_index[n5*64 + 9], arr_index[n5*64 + 11], n5%2==0);
        compareSwap(arr[n5*64 + 12], arr[n5*64 + 14], arr_index[n5*64 + 12], arr_index[n5*64 + 14], n5%2==0);
        compareSwap(arr[n5*64 + 13], arr[n5*64 + 15], arr_index[n5*64 + 13], arr_index[n5*64 + 15], n5%2==0);

        compareSwap(arr[n5*64 + 16], arr[n5*64 + 18], arr_index[n5*64 + 16], arr_index[n5*64 + 18], n5%2==0);
        compareSwap(arr[n5*64 + 17], arr[n5*64 + 19], arr_index[n5*64 + 17], arr_index[n5*64 + 19], n5%2==0);
        compareSwap(arr[n5*64 + 20], arr[n5*64 + 22], arr_index[n5*64 + 20], arr_index[n5*64 + 22], n5%2==0);
        compareSwap(arr[n5*64 + 21], arr[n5*64 + 23], arr_index[n5*64 + 21], arr_index[n5*64 + 23], n5%2==0);
        compareSwap(arr[n5*64 + 24], arr[n5*64 + 26], arr_index[n5*64 + 24], arr_index[n5*64 + 26], n5%2==0);
        compareSwap(arr[n5*64 + 25], arr[n5*64 + 27], arr_index[n5*64 + 25], arr_index[n5*64 + 27], n5%2==0);
        compareSwap(arr[n5*64 + 28], arr[n5*64 + 30], arr_index[n5*64 + 28], arr_index[n5*64 + 30], n5%2==0);
        compareSwap(arr[n5*64 + 29], arr[n5*64 + 31], arr_index[n5*64 + 29], arr_index[n5*64 + 31], n5%2==0);
        compareSwap(arr[n5*64 + 32], arr[n5*64 + 34], arr_index[n5*64 + 32], arr_index[n5*64 + 34], n5%2==0);
        compareSwap(arr[n5*64 + 33], arr[n5*64 + 35], arr_index[n5*64 + 33], arr_index[n5*64 + 35], n5%2==0);
        
        compareSwap(arr[n5*64 + 36], arr[n5*64 + 38], arr_index[n5*64 + 36], arr_index[n5*64 + 38], n5%2==0);
        compareSwap(arr[n5*64 + 37], arr[n5*64 + 39], arr_index[n5*64 + 37], arr_index[n5*64 + 39], n5%2==0);
        compareSwap(arr[n5*64 + 40], arr[n5*64 + 42], arr_index[n5*64 + 40], arr_index[n5*64 + 42], n5%2==0);
        compareSwap(arr[n5*64 + 41], arr[n5*64 + 43], arr_index[n5*64 + 41], arr_index[n5*64 + 43], n5%2==0);
        compareSwap(arr[n5*64 + 44], arr[n5*64 + 46], arr_index[n5*64 + 44], arr_index[n5*64 + 46], n5%2==0);
        compareSwap(arr[n5*64 + 45], arr[n5*64 + 47], arr_index[n5*64 + 45], arr_index[n5*64 + 47], n5%2==0);
        compareSwap(arr[n5*64 + 48], arr[n5*64 + 50], arr_index[n5*64 + 48], arr_index[n5*64 + 50], n5%2==0);
        compareSwap(arr[n5*64 + 49], arr[n5*64 + 51], arr_index[n5*64 + 49], arr_index[n5*64 + 51], n5%2==0);
        compareSwap(arr[n5*64 + 52], arr[n5*64 + 54], arr_index[n5*64 + 52], arr_index[n5*64 + 54], n5%2==0);
        compareSwap(arr[n5*64 + 53], arr[n5*64 + 55], arr_index[n5*64 + 53], arr_index[n5*64 + 55], n5%2==0);
        compareSwap(arr[n5*64 + 56], arr[n5*64 + 58], arr_index[n5*64 + 56], arr_index[n5*64 + 58], n5%2==0);
        compareSwap(arr[n5*64 + 57], arr[n5*64 + 59], arr_index[n5*64 + 57], arr_index[n5*64 + 59], n5%2==0);
        compareSwap(arr[n5*64 + 60], arr[n5*64 + 62], arr_index[n5*64 + 60], arr_index[n5*64 + 62], n5%2==0);
        compareSwap(arr[n5*64 + 61], arr[n5*64 + 63], arr_index[n5*64 + 61], arr_index[n5*64 + 63], n5%2==0);
    }

    for(int n6 = 0; n6 < 2; n6++){

        compareSwap(arr[n6*64], arr[n6*64 + 1], arr_index[n6*64], arr_index[n6*64 + 1], n6%2==0);
        compareSwap(arr[n6*64 + 2], arr[n6*64 + 3], arr_index[n6*64 + 2], arr_index[n6*64 + 3], n6%2==0);
        compareSwap(arr[n6*64 + 4], arr[n6*64 + 5], arr_index[n6*64 + 4], arr_index[n6*64 + 5], n6%2==0);
        compareSwap(arr[n6*64 + 6], arr[n6*64 + 7], arr_index[n6*64 + 6], arr_index[n6*64 + 7], n6%2==0);
        compareSwap(arr[n6*64 + 8], arr[n6*64 + 9], arr_index[n6*64 + 8], arr_index[n6*64 + 9], n6%2==0);
        compareSwap(arr[n6*64 + 10], arr[n6*64 + 11], arr_index[n6*64 + 10], arr_index[n6*64 + 11], n6%2==0);
        compareSwap(arr[n6*64 + 12], arr[n6*64 + 13], arr_index[n6*64 + 12], arr_index[n6*64 + 13], n6%2==0);
        compareSwap(arr[n6*64 + 14], arr[n6*64 + 15], arr_index[n6*64 + 14], arr_index[n6*64 + 15], n6%2==0);
        compareSwap(arr[n6*64 + 16], arr[n6*64 + 17], arr_index[n6*64 + 16], arr_index[n6*64 + 17], n6%2==0);
        compareSwap(arr[n6*64 + 18], arr[n6*64 + 19], arr_index[n6*64 + 18], arr_index[n6*64 + 19], n6%2==0);
        compareSwap(arr[n6*64 + 20], arr[n6*64 + 21], arr_index[n6*64 + 20], arr_index[n6*64 + 21], n6%2==0);
        compareSwap(arr[n6*64 + 22], arr[n6*64 + 23], arr_index[n6*64 + 22], arr_index[n6*64 + 23], n6%2==0);
        compareSwap(arr[n6*64 + 24], arr[n6*64 + 25], arr_index[n6*64 + 24], arr_index[n6*64 + 25], n6%2==0);
        compareSwap(arr[n6*64 + 26], arr[n6*64 + 27], arr_index[n6*64 + 26], arr_index[n6*64 + 27], n6%2==0);
        compareSwap(arr[n6*64 + 28], arr[n6*64 + 29], arr_index[n6*64 + 28], arr_index[n6*64 + 29], n6%2==0);
        compareSwap(arr[n6*64 + 30], arr[n6*64 + 31], arr_index[n6*64 + 30], arr_index[n6*64 + 31], n6%2==0);
        compareSwap(arr[n6*64 + 32], arr[n6*64 + 33], arr_index[n6*64 + 32], arr_index[n6*64 + 33], n6%2==0);
        compareSwap(arr[n6*64 + 34], arr[n6*64 + 35], arr_index[n6*64 + 34], arr_index[n6*64 + 35], n6%2==0);
        compareSwap(arr[n6*64 + 36], arr[n6*64 + 37], arr_index[n6*64 + 36], arr_index[n6*64 + 37], n6%2==0);
        compareSwap(arr[n6*64 + 38], arr[n6*64 + 39], arr_index[n6*64 + 38], arr_index[n6*64 + 39], n6%2==0);
        compareSwap(arr[n6*64 + 40], arr[n6*64 + 41], arr_index[n6*64 + 40], arr_index[n6*64 + 41], n6%2==0);
        compareSwap(arr[n6*64 + 42], arr[n6*64 + 43], arr_index[n6*64 + 42], arr_index[n6*64 + 43], n6%2==0);
        compareSwap(arr[n6*64 + 44], arr[n6*64 + 45], arr_index[n6*64 + 44], arr_index[n6*64 + 45], n6%2==0);
        compareSwap(arr[n6*64 + 46], arr[n6*64 + 47], arr_index[n6*64 + 46], arr_index[n6*64 + 47], n6%2==0);
        compareSwap(arr[n6*64 + 48], arr[n6*64 + 49], arr_index[n6*64 + 48], arr_index[n6*64 + 49], n6%2==0);
        compareSwap(arr[n6*64 + 50], arr[n6*64 + 51], arr_index[n6*64 + 50], arr_index[n6*64 + 51], n6%2==0);
        compareSwap(arr[n6*64 + 52], arr[n6*64 + 53], arr_index[n6*64 + 52], arr_index[n6*64 + 53], n6%2==0);
        compareSwap(arr[n6*64 + 54], arr[n6*64 + 55], arr_index[n6*64 + 54], arr_index[n6*64 + 55], n6%2==0);
        compareSwap(arr[n6*64 + 56], arr[n6*64 + 57], arr_index[n6*64 + 56], arr_index[n6*64 + 57], n6%2==0);
        compareSwap(arr[n6*64 + 58], arr[n6*64 + 59], arr_index[n6*64 + 58], arr_index[n6*64 + 59], n6%2==0);
        compareSwap(arr[n6*64 + 60], arr[n6*64 + 61], arr_index[n6*64 + 60], arr_index[n6*64 + 61], n6%2==0);
        compareSwap(arr[n6*64 + 62], arr[n6*64 + 63], arr_index[n6*64 + 62], arr_index[n6*64 + 63], n6%2==0);
    }
    //Stage 7
    for(int o = 0; o < 1; o++){

        compareSwap(arr[o*128], arr[o*128 + 64], arr_index[o*128], arr_index[o*128 + 64], o%2==0);
        compareSwap(arr[o*128 + 1], arr[o*128 + 65], arr_index[o*128 + 1], arr_index[o*128 + 65], o%2==0);
        compareSwap(arr[o*128 + 2], arr[o*128 + 66], arr_index[o*128 + 2], arr_index[o*128 + 66], o%2==0);
        compareSwap(arr[o*128 + 3], arr[o*128 + 67], arr_index[o*128 + 3], arr_index[o*128 + 67], o%2==0);
        compareSwap(arr[o*128 + 4], arr[o*128 + 68], arr_index[o*128 + 4], arr_index[o*128 + 68], o%2==0);
        compareSwap(arr[o*128 + 5], arr[o*128 + 69], arr_index[o*128 + 5], arr_index[o*128 + 69], o%2==0);
        compareSwap(arr[o*128 + 6], arr[o*128 + 70], arr_index[o*128 + 6], arr_index[o*128 + 70], o%2==0);
        compareSwap(arr[o*128 + 7], arr[o*128 + 71], arr_index[o*128 + 7], arr_index[o*128 + 71], o%2==0);
        compareSwap(arr[o*128 + 8], arr[o*128 + 72], arr_index[o*128 + 8], arr_index[o*128 + 72], o%2==0);
        compareSwap(arr[o*128 + 9], arr[o*128 + 73], arr_index[o*128 + 9], arr_index[o*128 + 73], o%2==0);
        compareSwap(arr[o*128 + 10], arr[o*128 + 74], arr_index[o*128 + 10], arr_index[o*128 + 74], o%2==0);
        compareSwap(arr[o*128 + 11], arr[o*128 + 75], arr_index[o*128 + 11], arr_index[o*128 + 75], o%2==0);
        compareSwap(arr[o*128 + 12], arr[o*128 + 76], arr_index[o*128 + 12], arr_index[o*128 + 76], o%2==0);
        compareSwap(arr[o*128 + 13], arr[o*128 + 77], arr_index[o*128 + 13], arr_index[o*128 + 77], o%2==0);
        compareSwap(arr[o*128 + 14], arr[o*128 + 78], arr_index[o*128 + 14], arr_index[o*128 + 78], o%2==0);
        compareSwap(arr[o*128 + 15], arr[o*128 + 79], arr_index[o*128 + 15], arr_index[o*128 + 79], o%2==0);
        compareSwap(arr[o*128 + 16], arr[o*128 + 80], arr_index[o*128 + 16], arr_index[o*128 + 80], o%2==0);
        compareSwap(arr[o*128 + 17], arr[o*128 + 81], arr_index[o*128 + 17], arr_index[o*128 + 81], o%2==0);
        compareSwap(arr[o*128 + 18], arr[o*128 + 82], arr_index[o*128 + 18], arr_index[o*128 + 82], o%2==0);
        compareSwap(arr[o*128 + 19], arr[o*128 + 83], arr_index[o*128 + 19], arr_index[o*128 + 83], o%2==0);
        compareSwap(arr[o*128 + 20], arr[o*128 + 84], arr_index[o*128 + 20], arr_index[o*128 + 84], o%2==0);
        compareSwap(arr[o*128 + 21], arr[o*128 + 85], arr_index[o*128 + 21], arr_index[o*128 + 85], o%2==0);
        compareSwap(arr[o*128 + 22], arr[o*128 + 86], arr_index[o*128 + 22], arr_index[o*128 + 86], o%2==0);
        compareSwap(arr[o*128 + 23], arr[o*128 + 87], arr_index[o*128 + 23], arr_index[o*128 + 87], o%2==0);
        compareSwap(arr[o*128 + 24], arr[o*128 + 88], arr_index[o*128 + 24], arr_index[o*128 + 88], o%2==0);
        compareSwap(arr[o*128 + 25], arr[o*128 + 89], arr_index[o*128 + 25], arr_index[o*128 + 89], o%2==0);
        compareSwap(arr[o*128 + 26], arr[o*128 + 90], arr_index[o*128 + 26], arr_index[o*128 + 90], o%2==0);
        compareSwap(arr[o*128 + 27], arr[o*128 + 91], arr_index[o*128 + 27], arr_index[o*128 + 91], o%2==0);
        compareSwap(arr[o*128 + 28], arr[o*128 + 92], arr_index[o*128 + 28], arr_index[o*128 + 92], o%2==0);
        compareSwap(arr[o*128 + 29], arr[o*128 + 93], arr_index[o*128 + 29], arr_index[o*128 + 93], o%2==0);
        compareSwap(arr[o*128 + 30], arr[o*128 + 94], arr_index[o*128 + 30], arr_index[o*128 + 94], o%2==0);
        compareSwap(arr[o*128 + 31], arr[o*128 + 95], arr_index[o*128 + 31], arr_index[o*128 + 95], o%2==0);
        compareSwap(arr[o*128 + 32], arr[o*128 + 96], arr_index[o*128 + 32], arr_index[o*128 + 96], o%2==0);
        compareSwap(arr[o*128 + 33], arr[o*128 + 97], arr_index[o*128 + 33], arr_index[o*128 + 97], o%2==0);
        compareSwap(arr[o*128 + 34], arr[o*128 + 98], arr_index[o*128 + 34], arr_index[o*128 + 98], o%2==0);
        compareSwap(arr[o*128 + 35], arr[o*128 + 99], arr_index[o*128 + 35], arr_index[o*128 + 99], o%2==0);
        compareSwap(arr[o*128 + 36], arr[o*128 + 100], arr_index[o*128 + 36], arr_index[o*128 + 100], o%2==0);
        compareSwap(arr[o*128 + 37], arr[o*128 + 101], arr_index[o*128 + 37], arr_index[o*128 + 101], o%2==0);
        compareSwap(arr[o*128 + 38], arr[o*128 + 102], arr_index[o*128 + 38], arr_index[o*128 + 102], o%2==0);
        compareSwap(arr[o*128 + 39], arr[o*128 + 103], arr_index[o*128 + 39], arr_index[o*128 + 103], o%2==0);
        compareSwap(arr[o*128 + 40], arr[o*128 + 104], arr_index[o*128 + 40], arr_index[o*128 + 104], o%2==0);
        compareSwap(arr[o*128 + 41], arr[o*128 + 105], arr_index[o*128 + 41], arr_index[o*128 + 105], o%2==0);
        compareSwap(arr[o*128 + 42], arr[o*128 + 106], arr_index[o*128 + 42], arr_index[o*128 + 106], o%2==0);
        compareSwap(arr[o*128 + 43], arr[o*128 + 107], arr_index[o*128 + 43], arr_index[o*128 + 107], o%2==0);
        compareSwap(arr[o*128 + 44], arr[o*128 + 108], arr_index[o*128 + 44], arr_index[o*128 + 108], o%2==0);
        compareSwap(arr[o*128 + 45], arr[o*128 + 109], arr_index[o*128 + 45], arr_index[o*128 + 109], o%2==0);
        compareSwap(arr[o*128 + 46], arr[o*128 + 110], arr_index[o*128 + 46], arr_index[o*128 + 110], o%2==0);
        compareSwap(arr[o*128 + 47], arr[o*128 + 111], arr_index[o*128 + 47], arr_index[o*128 + 111], o%2==0);
        compareSwap(arr[o*128 + 48], arr[o*128 + 112], arr_index[o*128 + 48], arr_index[o*128 + 112], o%2==0);
        compareSwap(arr[o*128 + 49], arr[o*128 + 113], arr_index[o*128 + 49], arr_index[o*128 + 113], o%2==0);
        compareSwap(arr[o*128 + 50], arr[o*128 + 114], arr_index[o*128 + 50], arr_index[o*128 + 114], o%2==0);
        compareSwap(arr[o*128 + 51], arr[o*128 + 115], arr_index[o*128 + 51], arr_index[o*128 + 115], o%2==0);
        compareSwap(arr[o*128 + 52], arr[o*128 + 116], arr_index[o*128 + 52], arr_index[o*128 + 116], o%2==0);
        compareSwap(arr[o*128 + 53], arr[o*128 + 117], arr_index[o*128 + 53], arr_index[o*128 + 117], o%2==0);
        compareSwap(arr[o*128 + 54], arr[o*128 + 118], arr_index[o*128 + 54], arr_index[o*128 + 118], o%2==0);
        compareSwap(arr[o*128 + 55], arr[o*128 + 119], arr_index[o*128 + 55], arr_index[o*128 + 119], o%2==0);
        compareSwap(arr[o*128 + 56], arr[o*128 + 120], arr_index[o*128 + 56], arr_index[o*128 + 120], o%2==0);
        compareSwap(arr[o*128 + 57], arr[o*128 + 121], arr_index[o*128 + 57], arr_index[o*128 + 121], o%2==0);
        compareSwap(arr[o*128 + 58], arr[o*128 + 122], arr_index[o*128 + 58], arr_index[o*128 + 122], o%2==0);
        compareSwap(arr[o*128 + 59], arr[o*128 + 123], arr_index[o*128 + 59], arr_index[o*128 + 123], o%2==0);
        compareSwap(arr[o*128 + 60], arr[o*128 + 124], arr_index[o*128 + 60], arr_index[o*128 + 124], o%2==0);
        compareSwap(arr[o*128 + 61], arr[o*128 + 125], arr_index[o*128 + 61], arr_index[o*128 + 125], o%2==0);
        compareSwap(arr[o*128 + 62], arr[o*128 + 126], arr_index[o*128 + 62], arr_index[o*128 + 126], o%2==0);
        compareSwap(arr[o*128 + 63], arr[o*128 + 127], arr_index[o*128 + 63], arr_index[o*128 + 127], o%2==0);
    }
    for(int o1 = 0; o1 < 1; o1++){

        compareSwap(arr[o1*128], arr[o1*128 + 32], arr_index[o1*128], arr_index[o1*128 + 32], o1%2==0);
        compareSwap(arr[o1*128 + 1], arr[o1*128 + 33], arr_index[o1*128 + 1], arr_index[o1*128 + 33], o1%2==0);
        compareSwap(arr[o1*128 + 2], arr[o1*128 + 34], arr_index[o1*128 + 2], arr_index[o1*128 + 34], o1%2==0);
        compareSwap(arr[o1*128 + 3], arr[o1*128 + 35], arr_index[o1*128 + 3], arr_index[o1*128 + 35], o1%2==0);
        compareSwap(arr[o1*128 + 4], arr[o1*128 + 36], arr_index[o1*128 + 4], arr_index[o1*128 + 36], o1%2==0);
        compareSwap(arr[o1*128 + 5], arr[o1*128 + 37], arr_index[o1*128 + 5], arr_index[o1*128 + 37], o1%2==0);
        compareSwap(arr[o1*128 + 6], arr[o1*128 + 38], arr_index[o1*128 + 6], arr_index[o1*128 + 38], o1%2==0);
        compareSwap(arr[o1*128 + 7], arr[o1*128 + 39], arr_index[o1*128 + 7], arr_index[o1*128 + 39], o1%2==0);
        compareSwap(arr[o1*128 + 8], arr[o1*128 + 40], arr_index[o1*128 + 8], arr_index[o1*128 + 40], o1%2==0);
        compareSwap(arr[o1*128 + 9], arr[o1*128 + 41], arr_index[o1*128 + 9], arr_index[o1*128 + 41], o1%2==0);
        compareSwap(arr[o1*128 + 10], arr[o1*128 + 42], arr_index[o1*128 + 10], arr_index[o1*128 + 42], o1%2==0);
        compareSwap(arr[o1*128 + 11], arr[o1*128 + 43], arr_index[o1*128 + 11], arr_index[o1*128 + 43], o1%2==0);
        compareSwap(arr[o1*128 + 12], arr[o1*128 + 44], arr_index[o1*128 + 12], arr_index[o1*128 + 44], o1%2==0);
        compareSwap(arr[o1*128 + 13], arr[o1*128 + 45], arr_index[o1*128 + 13], arr_index[o1*128 + 45], o1%2==0);
        compareSwap(arr[o1*128 + 14], arr[o1*128 + 46], arr_index[o1*128 + 14], arr_index[o1*128 + 46], o1%2==0);
        compareSwap(arr[o1*128 + 15], arr[o1*128 + 47], arr_index[o1*128 + 15], arr_index[o1*128 + 47], o1%2==0);
        compareSwap(arr[o1*128 + 16], arr[o1*128 + 48], arr_index[o1*128 + 16], arr_index[o1*128 + 48], o1%2==0);
        compareSwap(arr[o1*128 + 17], arr[o1*128 + 49], arr_index[o1*128 + 17], arr_index[o1*128 + 49], o1%2==0);
        compareSwap(arr[o1*128 + 18], arr[o1*128 + 50], arr_index[o1*128 + 18], arr_index[o1*128 + 50], o1%2==0);
        compareSwap(arr[o1*128 + 19], arr[o1*128 + 51], arr_index[o1*128 + 19], arr_index[o1*128 + 51], o1%2==0);
        compareSwap(arr[o1*128 + 20], arr[o1*128 + 52], arr_index[o1*128 + 20], arr_index[o1*128 + 52], o1%2==0);
        compareSwap(arr[o1*128 + 21], arr[o1*128 + 53], arr_index[o1*128 + 21], arr_index[o1*128 + 53], o1%2==0);
        compareSwap(arr[o1*128 + 22], arr[o1*128 + 54], arr_index[o1*128 + 22], arr_index[o1*128 + 54], o1%2==0);
        compareSwap(arr[o1*128 + 23], arr[o1*128 + 55], arr_index[o1*128 + 23], arr_index[o1*128 + 55], o1%2==0);
        compareSwap(arr[o1*128 + 24], arr[o1*128 + 56], arr_index[o1*128 + 24], arr_index[o1*128 + 56], o1%2==0);
        compareSwap(arr[o1*128 + 25], arr[o1*128 + 57], arr_index[o1*128 + 25], arr_index[o1*128 + 57], o1%2==0);
        compareSwap(arr[o1*128 + 26], arr[o1*128 + 58], arr_index[o1*128 + 26], arr_index[o1*128 + 58], o1%2==0);
        compareSwap(arr[o1*128 + 27], arr[o1*128 + 59], arr_index[o1*128 + 27], arr_index[o1*128 + 59], o1%2==0);
        compareSwap(arr[o1*128 + 28], arr[o1*128 + 60], arr_index[o1*128 + 28], arr_index[o1*128 + 60], o1%2==0);
        compareSwap(arr[o1*128 + 29], arr[o1*128 + 61], arr_index[o1*128 + 29], arr_index[o1*128 + 61], o1%2==0);
        compareSwap(arr[o1*128 + 30], arr[o1*128 + 62], arr_index[o1*128 + 30], arr_index[o1*128 + 62], o1%2==0);
        compareSwap(arr[o1*128 + 31], arr[o1*128 + 63], arr_index[o1*128 + 31], arr_index[o1*128 + 63], o1%2==0);
        
        compareSwap(arr[o1*128 + 64], arr[o1*128 + 96], arr_index[o1*128 + 64], arr_index[o1*128 + 96], o1%2==0);
        compareSwap(arr[o1*128 + 65], arr[o1*128 + 97], arr_index[o1*128 + 65], arr_index[o1*128 + 97], o1%2==0);
        compareSwap(arr[o1*128 + 66], arr[o1*128 + 98], arr_index[o1*128 + 66], arr_index[o1*128 + 98], o1%2==0);
        compareSwap(arr[o1*128 + 67], arr[o1*128 + 99], arr_index[o1*128 + 67], arr_index[o1*128 + 99], o1%2==0);
        compareSwap(arr[o1*128 + 68], arr[o1*128 + 100], arr_index[o1*128 + 68], arr_index[o1*128 + 100], o1%2==0);
        compareSwap(arr[o1*128 + 69], arr[o1*128 + 101], arr_index[o1*128 + 69], arr_index[o1*128 + 101], o1%2==0);
        compareSwap(arr[o1*128 + 70], arr[o1*128 + 102], arr_index[o1*128 + 70], arr_index[o1*128 + 102], o1%2==0);
        compareSwap(arr[o1*128 + 71], arr[o1*128 + 103], arr_index[o1*128 + 71], arr_index[o1*128 + 103], o1%2==0);
        compareSwap(arr[o1*128 + 72], arr[o1*128 + 104], arr_index[o1*128 + 72], arr_index[o1*128 + 104], o1%2==0);
        compareSwap(arr[o1*128 + 73], arr[o1*128 + 105], arr_index[o1*128 + 73], arr_index[o1*128 + 105], o1%2==0);
        compareSwap(arr[o1*128 + 74], arr[o1*128 + 106], arr_index[o1*128 + 74], arr_index[o1*128 + 106], o1%2==0);
        compareSwap(arr[o1*128 + 75], arr[o1*128 + 107], arr_index[o1*128 + 75], arr_index[o1*128 + 107], o1%2==0);
        compareSwap(arr[o1*128 + 76], arr[o1*128 + 108], arr_index[o1*128 + 76], arr_index[o1*128 + 108], o1%2==0);
        compareSwap(arr[o1*128 + 77], arr[o1*128 + 109], arr_index[o1*128 + 77], arr_index[o1*128 + 109], o1%2==0);
        compareSwap(arr[o1*128 + 78], arr[o1*128 + 110], arr_index[o1*128 + 78], arr_index[o1*128 + 110], o1%2==0);
        compareSwap(arr[o1*128 + 79], arr[o1*128 + 111], arr_index[o1*128 + 79], arr_index[o1*128 + 111], o1%2==0);
        compareSwap(arr[o1*128 + 80], arr[o1*128 + 112], arr_index[o1*128 + 80], arr_index[o1*128 + 112], o1%2==0);
        compareSwap(arr[o1*128 + 81], arr[o1*128 + 113], arr_index[o1*128 + 81], arr_index[o1*128 + 113], o1%2==0);
        compareSwap(arr[o1*128 + 82], arr[o1*128 + 114], arr_index[o1*128 + 82], arr_index[o1*128 + 114], o1%2==0);
        compareSwap(arr[o1*128 + 83], arr[o1*128 + 115], arr_index[o1*128 + 83], arr_index[o1*128 + 115], o1%2==0);
        compareSwap(arr[o1*128 + 84], arr[o1*128 + 116], arr_index[o1*128 + 84], arr_index[o1*128 + 116], o1%2==0);
        compareSwap(arr[o1*128 + 85], arr[o1*128 + 117], arr_index[o1*128 + 85], arr_index[o1*128 + 117], o1%2==0);
        compareSwap(arr[o1*128 + 86], arr[o1*128 + 118], arr_index[o1*128 + 86], arr_index[o1*128 + 118], o1%2==0);
        compareSwap(arr[o1*128 + 87], arr[o1*128 + 119], arr_index[o1*128 + 87], arr_index[o1*128 + 119], o1%2==0);
        compareSwap(arr[o1*128 + 88], arr[o1*128 + 120], arr_index[o1*128 + 88], arr_index[o1*128 + 120], o1%2==0);
        compareSwap(arr[o1*128 + 89], arr[o1*128 + 121], arr_index[o1*128 + 89], arr_index[o1*128 + 121], o1%2==0);
        compareSwap(arr[o1*128 + 90], arr[o1*128 + 122], arr_index[o1*128 + 90], arr_index[o1*128 + 122], o1%2==0);
        compareSwap(arr[o1*128 + 91], arr[o1*128 + 123], arr_index[o1*128 + 91], arr_index[o1*128 + 123], o1%2==0);
        compareSwap(arr[o1*128 + 92], arr[o1*128 + 124], arr_index[o1*128 + 92], arr_index[o1*128 + 124], o1%2==0);
        compareSwap(arr[o1*128 + 93], arr[o1*128 + 125], arr_index[o1*128 + 93], arr_index[o1*128 + 125], o1%2==0);
        compareSwap(arr[o1*128 + 94], arr[o1*128 + 126], arr_index[o1*128 + 94], arr_index[o1*128 + 126], o1%2==0);
        compareSwap(arr[o1*128 + 95], arr[o1*128 + 127], arr_index[o1*128 + 95], arr_index[o1*128 + 127], o1%2==0);
    }
    for(int o2 = 0; o2 < 1; o2++){

        compareSwap(arr[o2*128], arr[o2*128 + 16], arr_index[o2*128], arr_index[o2*128 + 16], o2%2==0);
        compareSwap(arr[o2*128 + 1], arr[o2*128 + 17], arr_index[o2*128 + 1], arr_index[o2*128 + 17], o2%2==0);
        compareSwap(arr[o2*128 + 2], arr[o2*128 + 18], arr_index[o2*128 + 2], arr_index[o2*128 + 18], o2%2==0);
        compareSwap(arr[o2*128 + 3], arr[o2*128 + 19], arr_index[o2*128 + 3], arr_index[o2*128 + 19], o2%2==0);
        compareSwap(arr[o2*128 + 4], arr[o2*128 + 20], arr_index[o2*128 + 4], arr_index[o2*128 + 20], o2%2==0);
        compareSwap(arr[o2*128 + 5], arr[o2*128 + 21], arr_index[o2*128 + 5], arr_index[o2*128 + 21], o2%2==0);
        compareSwap(arr[o2*128 + 6], arr[o2*128 + 22], arr_index[o2*128 + 6], arr_index[o2*128 + 22], o2%2==0);
        compareSwap(arr[o2*128 + 7], arr[o2*128 + 23], arr_index[o2*128 + 7], arr_index[o2*128 + 23], o2%2==0);
        compareSwap(arr[o2*128 + 8], arr[o2*128 + 24], arr_index[o2*128 + 8], arr_index[o2*128 + 24], o2%2==0);
        compareSwap(arr[o2*128 + 9], arr[o2*128 + 25], arr_index[o2*128 + 9], arr_index[o2*128 + 25], o2%2==0);
        compareSwap(arr[o2*128 + 10], arr[o2*128 + 26], arr_index[o2*128 + 10], arr_index[o2*128 + 26], o2%2==0);
        compareSwap(arr[o2*128 + 11], arr[o2*128 + 27], arr_index[o2*128 + 11], arr_index[o2*128 + 27], o2%2==0);
        compareSwap(arr[o2*128 + 12], arr[o2*128 + 28], arr_index[o2*128 + 12], arr_index[o2*128 + 28], o2%2==0);
        compareSwap(arr[o2*128 + 13], arr[o2*128 + 29], arr_index[o2*128 + 13], arr_index[o2*128 + 29], o2%2==0);
        compareSwap(arr[o2*128 + 14], arr[o2*128 + 30], arr_index[o2*128 + 14], arr_index[o2*128 + 30], o2%2==0);
        compareSwap(arr[o2*128 + 15], arr[o2*128 + 31], arr_index[o2*128 + 15], arr_index[o2*128 + 31], o2%2==0);

        compareSwap(arr[o2*128 + 32], arr[o2*128 + 48], arr_index[o2*128 + 32], arr_index[o2*128 + 48], o2%2==0);
        compareSwap(arr[o2*128 + 33], arr[o2*128 + 49], arr_index[o2*128 + 33], arr_index[o2*128 + 49], o2%2==0);
        compareSwap(arr[o2*128 + 34], arr[o2*128 + 50], arr_index[o2*128 + 34], arr_index[o2*128 + 50], o2%2==0);
        compareSwap(arr[o2*128 + 35], arr[o2*128 + 51], arr_index[o2*128 + 35], arr_index[o2*128 + 51], o2%2==0);
        compareSwap(arr[o2*128 + 36], arr[o2*128 + 52], arr_index[o2*128 + 36], arr_index[o2*128 + 52], o2%2==0);
        compareSwap(arr[o2*128 + 37], arr[o2*128 + 53], arr_index[o2*128 + 37], arr_index[o2*128 + 53], o2%2==0);
        compareSwap(arr[o2*128 + 38], arr[o2*128 + 54], arr_index[o2*128 + 38], arr_index[o2*128 + 54], o2%2==0);
        compareSwap(arr[o2*128 + 39], arr[o2*128 + 55], arr_index[o2*128 + 39], arr_index[o2*128 + 55], o2%2==0);
        compareSwap(arr[o2*128 + 40], arr[o2*128 + 56], arr_index[o2*128 + 40], arr_index[o2*128 + 56], o2%2==0);
        compareSwap(arr[o2*128 + 41], arr[o2*128 + 57], arr_index[o2*128 + 41], arr_index[o2*128 + 57], o2%2==0);
        compareSwap(arr[o2*128 + 42], arr[o2*128 + 58], arr_index[o2*128 + 42], arr_index[o2*128 + 58], o2%2==0);
        compareSwap(arr[o2*128 + 43], arr[o2*128 + 59], arr_index[o2*128 + 43], arr_index[o2*128 + 59], o2%2==0);
        compareSwap(arr[o2*128 + 44], arr[o2*128 + 60], arr_index[o2*128 + 44], arr_index[o2*128 + 60], o2%2==0);
        compareSwap(arr[o2*128 + 45], arr[o2*128 + 61], arr_index[o2*128 + 45], arr_index[o2*128 + 61], o2%2==0);
        compareSwap(arr[o2*128 + 46], arr[o2*128 + 62], arr_index[o2*128 + 46], arr_index[o2*128 + 62], o2%2==0);
        compareSwap(arr[o2*128 + 47], arr[o2*128 + 63], arr_index[o2*128 + 47], arr_index[o2*128 + 63], o2%2==0);

        compareSwap(arr[o2*128 + 64], arr[o2*128 + 80], arr_index[o2*128 + 64], arr_index[o2*128 + 80], o2%2==0);
        compareSwap(arr[o2*128 + 65], arr[o2*128 + 81], arr_index[o2*128 + 65], arr_index[o2*128 + 81], o2%2==0);
        compareSwap(arr[o2*128 + 66], arr[o2*128 + 82], arr_index[o2*128 + 66], arr_index[o2*128 + 82], o2%2==0);
        compareSwap(arr[o2*128 + 67], arr[o2*128 + 83], arr_index[o2*128 + 67], arr_index[o2*128 + 83], o2%2==0);
        compareSwap(arr[o2*128 + 68], arr[o2*128 + 84], arr_index[o2*128 + 68], arr_index[o2*128 + 84], o2%2==0);
        compareSwap(arr[o2*128 + 69], arr[o2*128 + 85], arr_index[o2*128 + 69], arr_index[o2*128 + 85], o2%2==0);
        compareSwap(arr[o2*128 + 70], arr[o2*128 + 86], arr_index[o2*128 + 70], arr_index[o2*128 + 86], o2%2==0);
        compareSwap(arr[o2*128 + 71], arr[o2*128 + 87], arr_index[o2*128 + 71], arr_index[o2*128 + 87], o2%2==0);
        compareSwap(arr[o2*128 + 72], arr[o2*128 + 88], arr_index[o2*128 + 72], arr_index[o2*128 + 88], o2%2==0);
        compareSwap(arr[o2*128 + 73], arr[o2*128 + 89], arr_index[o2*128 + 73], arr_index[o2*128 + 89], o2%2==0);
        compareSwap(arr[o2*128 + 74], arr[o2*128 + 90], arr_index[o2*128 + 74], arr_index[o2*128 + 90], o2%2==0);
        compareSwap(arr[o2*128 + 75], arr[o2*128 + 91], arr_index[o2*128 + 75], arr_index[o2*128 + 91], o2%2==0);
        compareSwap(arr[o2*128 + 76], arr[o2*128 + 92], arr_index[o2*128 + 76], arr_index[o2*128 + 92], o2%2==0);
        compareSwap(arr[o2*128 + 77], arr[o2*128 + 93], arr_index[o2*128 + 77], arr_index[o2*128 + 93], o2%2==0);
        compareSwap(arr[o2*128 + 78], arr[o2*128 + 94], arr_index[o2*128 + 78], arr_index[o2*128 + 94], o2%2==0);
        compareSwap(arr[o2*128 + 79], arr[o2*128 + 95], arr_index[o2*128 + 79], arr_index[o2*128 + 95], o2%2==0);

        compareSwap(arr[o2*128 + 96], arr[o2*128 + 112], arr_index[o2*128 + 96], arr_index[o2*128 + 112], o2%2==0);
        compareSwap(arr[o2*128 + 97], arr[o2*128 + 113], arr_index[o2*128 + 97], arr_index[o2*128 + 113], o2%2==0);
        compareSwap(arr[o2*128 + 98], arr[o2*128 + 114], arr_index[o2*128 + 98], arr_index[o2*128 + 114], o2%2==0);
        compareSwap(arr[o2*128 + 99], arr[o2*128 + 115], arr_index[o2*128 + 99], arr_index[o2*128 + 115], o2%2==0);
        compareSwap(arr[o2*128 + 100], arr[o2*128 + 116], arr_index[o2*128 + 100], arr_index[o2*128 + 116], o2%2==0);
        compareSwap(arr[o2*128 + 101], arr[o2*128 + 117], arr_index[o2*128 + 101], arr_index[o2*128 + 117], o2%2==0);
        compareSwap(arr[o2*128 + 102], arr[o2*128 + 118], arr_index[o2*128 + 102], arr_index[o2*128 + 118], o2%2==0);
        compareSwap(arr[o2*128 + 103], arr[o2*128 + 119], arr_index[o2*128 + 103], arr_index[o2*128 + 119], o2%2==0);
        compareSwap(arr[o2*128 + 104], arr[o2*128 + 120], arr_index[o2*128 + 104], arr_index[o2*128 + 120], o2%2==0);
        compareSwap(arr[o2*128 + 105], arr[o2*128 + 121], arr_index[o2*128 + 105], arr_index[o2*128 + 121], o2%2==0);
        compareSwap(arr[o2*128 + 106], arr[o2*128 + 122], arr_index[o2*128 + 106], arr_index[o2*128 + 122], o2%2==0);
        compareSwap(arr[o2*128 + 107], arr[o2*128 + 123], arr_index[o2*128 + 107], arr_index[o2*128 + 123], o2%2==0);
        compareSwap(arr[o2*128 + 108], arr[o2*128 + 124], arr_index[o2*128 + 108], arr_index[o2*128 + 124], o2%2==0);
        compareSwap(arr[o2*128 + 109], arr[o2*128 + 125], arr_index[o2*128 + 109], arr_index[o2*128 + 125], o2%2==0);
        compareSwap(arr[o2*128 + 110], arr[o2*128 + 126], arr_index[o2*128 + 110], arr_index[o2*128 + 126], o2%2==0);
        compareSwap(arr[o2*128 + 111], arr[o2*128 + 127], arr_index[o2*128 + 111], arr_index[o2*128 + 127], o2%2==0);
    }
    for(int o3 = 0; o3 < 1; o3++){

        compareSwap(arr[o3*128], arr[o3*128 + 8], arr_index[o3*128], arr_index[o3*128 + 8], o3%2==0);
        compareSwap(arr[o3*128 + 1], arr[o3*128 + 9], arr_index[o3*128 + 1], arr_index[o3*128 + 9], o3%2==0);
        compareSwap(arr[o3*128 + 2], arr[o3*128 + 10], arr_index[o3*128 + 2], arr_index[o3*128 + 10], o3%2==0);
        compareSwap(arr[o3*128 + 3], arr[o3*128 + 11], arr_index[o3*128 + 3], arr_index[o3*128 + 11], o3%2==0);
        compareSwap(arr[o3*128 + 4], arr[o3*128 + 12], arr_index[o3*128 + 4], arr_index[o3*128 + 12], o3%2==0);
        compareSwap(arr[o3*128 + 5], arr[o3*128 + 13], arr_index[o3*128 + 5], arr_index[o3*128 + 13], o3%2==0);
        compareSwap(arr[o3*128 + 6], arr[o3*128 + 14], arr_index[o3*128 + 6], arr_index[o3*128 + 14], o3%2==0);
        compareSwap(arr[o3*128 + 7], arr[o3*128 + 15], arr_index[o3*128 + 7], arr_index[o3*128 + 15], o3%2==0);

        compareSwap(arr[o3*128 + 16], arr[o3*128 + 24], arr_index[o3*128 + 16], arr_index[o3*128 + 24], o3%2==0);
        compareSwap(arr[o3*128 + 17], arr[o3*128 + 25], arr_index[o3*128 + 17], arr_index[o3*128 + 25], o3%2==0);
        compareSwap(arr[o3*128 + 18], arr[o3*128 + 26], arr_index[o3*128 + 18], arr_index[o3*128 + 26], o3%2==0);
        compareSwap(arr[o3*128 + 19], arr[o3*128 + 27], arr_index[o3*128 + 19], arr_index[o3*128 + 27], o3%2==0);
        compareSwap(arr[o3*128 + 20], arr[o3*128 + 28], arr_index[o3*128 + 20], arr_index[o3*128 + 28], o3%2==0);
        compareSwap(arr[o3*128 + 21], arr[o3*128 + 29], arr_index[o3*128 + 21], arr_index[o3*128 + 29], o3%2==0);
        compareSwap(arr[o3*128 + 22], arr[o3*128 + 30], arr_index[o3*128 + 22], arr_index[o3*128 + 30], o3%2==0);
        compareSwap(arr[o3*128 + 23], arr[o3*128 + 31], arr_index[o3*128 + 23], arr_index[o3*128 + 31], o3%2==0);

        compareSwap(arr[o3*128 + 32], arr[o3*128 + 40], arr_index[o3*128 + 32], arr_index[o3*128 + 40], o3%2==0);
        compareSwap(arr[o3*128 + 33], arr[o3*128 + 41], arr_index[o3*128 + 33], arr_index[o3*128 + 41], o3%2==0);
        compareSwap(arr[o3*128 + 34], arr[o3*128 + 42], arr_index[o3*128 + 34], arr_index[o3*128 + 42], o3%2==0);
        compareSwap(arr[o3*128 + 35], arr[o3*128 + 43], arr_index[o3*128 + 35], arr_index[o3*128 + 43], o3%2==0);
        compareSwap(arr[o3*128 + 36], arr[o3*128 + 44], arr_index[o3*128 + 36], arr_index[o3*128 + 44], o3%2==0);
        compareSwap(arr[o3*128 + 37], arr[o3*128 + 45], arr_index[o3*128 + 37], arr_index[o3*128 + 45], o3%2==0);
        compareSwap(arr[o3*128 + 38], arr[o3*128 + 46], arr_index[o3*128 + 38], arr_index[o3*128 + 46], o3%2==0);
        compareSwap(arr[o3*128 + 39], arr[o3*128 + 47], arr_index[o3*128 + 39], arr_index[o3*128 + 47], o3%2==0);

        compareSwap(arr[o3*128 + 48], arr[o3*128 + 56], arr_index[o3*128 + 48], arr_index[o3*128 + 56], o3%2==0);
        compareSwap(arr[o3*128 + 49], arr[o3*128 + 57], arr_index[o3*128 + 49], arr_index[o3*128 + 57], o3%2==0);
        compareSwap(arr[o3*128 + 50], arr[o3*128 + 58], arr_index[o3*128 + 50], arr_index[o3*128 + 58], o3%2==0);
        compareSwap(arr[o3*128 + 51], arr[o3*128 + 59], arr_index[o3*128 + 51], arr_index[o3*128 + 59], o3%2==0);
        compareSwap(arr[o3*128 + 52], arr[o3*128 + 60], arr_index[o3*128 + 52], arr_index[o3*128 + 60], o3%2==0);
        compareSwap(arr[o3*128 + 53], arr[o3*128 + 61], arr_index[o3*128 + 53], arr_index[o3*128 + 61], o3%2==0);
        compareSwap(arr[o3*128 + 54], arr[o3*128 + 62], arr_index[o3*128 + 54], arr_index[o3*128 + 62], o3%2==0);
        compareSwap(arr[o3*128 + 55], arr[o3*128 + 63], arr_index[o3*128 + 55], arr_index[o3*128 + 63], o3%2==0);

        compareSwap(arr[o3*128 + 64], arr[o3*128 + 72], arr_index[o3*128 + 64], arr_index[o3*128 + 72], o3%2==0);
        compareSwap(arr[o3*128 + 65], arr[o3*128 + 73], arr_index[o3*128 + 65], arr_index[o3*128 + 73], o3%2==0);
        compareSwap(arr[o3*128 + 66], arr[o3*128 + 74], arr_index[o3*128 + 66], arr_index[o3*128 + 74], o3%2==0);
        compareSwap(arr[o3*128 + 67], arr[o3*128 + 75], arr_index[o3*128 + 67], arr_index[o3*128 + 75], o3%2==0);
        compareSwap(arr[o3*128 + 68], arr[o3*128 + 76], arr_index[o3*128 + 68], arr_index[o3*128 + 76], o3%2==0);
        compareSwap(arr[o3*128 + 69], arr[o3*128 + 77], arr_index[o3*128 + 69], arr_index[o3*128 + 77], o3%2==0);
        compareSwap(arr[o3*128 + 70], arr[o3*128 + 78], arr_index[o3*128 + 70], arr_index[o3*128 + 78], o3%2==0);
        compareSwap(arr[o3*128 + 71], arr[o3*128 + 79], arr_index[o3*128 + 71], arr_index[o3*128 + 79], o3%2==0);

        compareSwap(arr[o3*128 + 80], arr[o3*128 + 88], arr_index[o3*128 + 80], arr_index[o3*128 + 88], o3%2==0);
        compareSwap(arr[o3*128 + 81], arr[o3*128 + 89], arr_index[o3*128 + 81], arr_index[o3*128 + 89], o3%2==0);
        compareSwap(arr[o3*128 + 82], arr[o3*128 + 90], arr_index[o3*128 + 82], arr_index[o3*128 + 90], o3%2==0);
        compareSwap(arr[o3*128 + 83], arr[o3*128 + 91], arr_index[o3*128 + 83], arr_index[o3*128 + 91], o3%2==0);
        compareSwap(arr[o3*128 + 84], arr[o3*128 + 92], arr_index[o3*128 + 84], arr_index[o3*128 + 92], o3%2==0);
        compareSwap(arr[o3*128 + 85], arr[o3*128 + 93], arr_index[o3*128 + 85], arr_index[o3*128 + 93], o3%2==0);
        compareSwap(arr[o3*128 + 86], arr[o3*128 + 94], arr_index[o3*128 + 86], arr_index[o3*128 + 94], o3%2==0);
        compareSwap(arr[o3*128 + 87], arr[o3*128 + 95], arr_index[o3*128 + 87], arr_index[o3*128 + 95], o3%2==0);

        compareSwap(arr[o3*128 + 96], arr[o3*128 + 104], arr_index[o3*128 + 96], arr_index[o3*128 + 104], o3%2==0);
        compareSwap(arr[o3*128 + 97], arr[o3*128 + 105], arr_index[o3*128 + 97], arr_index[o3*128 + 105], o3%2==0);
        compareSwap(arr[o3*128 + 98], arr[o3*128 + 106], arr_index[o3*128 + 98], arr_index[o3*128 + 106], o3%2==0);
        compareSwap(arr[o3*128 + 99], arr[o3*128 + 107], arr_index[o3*128 + 99], arr_index[o3*128 + 107], o3%2==0);
        compareSwap(arr[o3*128 + 100], arr[o3*128 + 108], arr_index[o3*128 + 100], arr_index[o3*128 + 108], o3%2==0);
        compareSwap(arr[o3*128 + 101], arr[o3*128 + 109], arr_index[o3*128 + 101], arr_index[o3*128 + 109], o3%2==0);
        compareSwap(arr[o3*128 + 102], arr[o3*128 + 110], arr_index[o3*128 + 102], arr_index[o3*128 + 110], o3%2==0);
        compareSwap(arr[o3*128 + 103], arr[o3*128 + 111], arr_index[o3*128 + 103], arr_index[o3*128 + 111], o3%2==0);

        compareSwap(arr[o3*128 + 112], arr[o3*128 + 120], arr_index[o3*128 + 112], arr_index[o3*128 + 120], o3%2==0);
        compareSwap(arr[o3*128 + 113], arr[o3*128 + 121], arr_index[o3*128 + 113], arr_index[o3*128 + 121], o3%2==0);
        compareSwap(arr[o3*128 + 114], arr[o3*128 + 122], arr_index[o3*128 + 114], arr_index[o3*128 + 122], o3%2==0);
        compareSwap(arr[o3*128 + 115], arr[o3*128 + 123], arr_index[o3*128 + 115], arr_index[o3*128 + 123], o3%2==0);
        compareSwap(arr[o3*128 + 116], arr[o3*128 + 124], arr_index[o3*128 + 116], arr_index[o3*128 + 124], o3%2==0);
        compareSwap(arr[o3*128 + 117], arr[o3*128 + 125], arr_index[o3*128 + 117], arr_index[o3*128 + 125], o3%2==0);
        compareSwap(arr[o3*128 + 118], arr[o3*128 + 126], arr_index[o3*128 + 118], arr_index[o3*128 + 126], o3%2==0);
        compareSwap(arr[o3*128 + 119], arr[o3*128 + 127], arr_index[o3*128 + 119], arr_index[o3*128 + 127], o3%2==0);
    }
    for(int o4 = 0; o4 < 1; o4++){

        compareSwap(arr[o4*128], arr[o4*128 + 4], arr_index[o4*128], arr_index[o4*128 + 4], o4%2==0);
        compareSwap(arr[o4*128 + 1], arr[o4*128 + 5], arr_index[o4*128 + 1], arr_index[o4*128 + 5], o4%2==0);
        compareSwap(arr[o4*128 + 2], arr[o4*128 + 6], arr_index[o4*128 + 2], arr_index[o4*128 + 6], o4%2==0);
        compareSwap(arr[o4*128 + 3], arr[o4*128 + 7], arr_index[o4*128 + 3], arr_index[o4*128 + 7], o4%2==0);

        compareSwap(arr[o4*128 + 8], arr[o4*128 + 12], arr_index[o4*128 + 8], arr_index[o4*128 + 12], o4%2==0);
        compareSwap(arr[o4*128 + 9], arr[o4*128 + 13], arr_index[o4*128 + 9], arr_index[o4*128 + 13], o4%2==0);
        compareSwap(arr[o4*128 + 10], arr[o4*128 + 14], arr_index[o4*128 + 10], arr_index[o4*128 + 14], o4%2==0);
        compareSwap(arr[o4*128 + 11], arr[o4*128 + 15], arr_index[o4*128 + 11], arr_index[o4*128 + 15], o4%2==0);

        compareSwap(arr[o4*128 + 16], arr[o4*128 + 20], arr_index[o4*128 + 16], arr_index[o4*128 + 20], o4%2==0);
        compareSwap(arr[o4*128 + 17], arr[o4*128 + 21], arr_index[o4*128 + 17], arr_index[o4*128 + 21], o4%2==0);
        compareSwap(arr[o4*128 + 18], arr[o4*128 + 22], arr_index[o4*128 + 18], arr_index[o4*128 + 22], o4%2==0);
        compareSwap(arr[o4*128 + 19], arr[o4*128 + 23], arr_index[o4*128 + 19], arr_index[o4*128 + 23], o4%2==0);

        compareSwap(arr[o4*128 + 24], arr[o4*128 + 28], arr_index[o4*128 + 24], arr_index[o4*128 + 28], o4%2==0);
        compareSwap(arr[o4*128 + 25], arr[o4*128 + 29], arr_index[o4*128 + 25], arr_index[o4*128 + 29], o4%2==0);
        compareSwap(arr[o4*128 + 26], arr[o4*128 + 30], arr_index[o4*128 + 26], arr_index[o4*128 + 30], o4%2==0);
        compareSwap(arr[o4*128 + 27], arr[o4*128 + 31], arr_index[o4*128 + 27], arr_index[o4*128 + 31], o4%2==0);

        compareSwap(arr[o4*128 + 32], arr[o4*128 + 36], arr_index[o4*128 + 32], arr_index[o4*128 + 36], o4%2==0);
        compareSwap(arr[o4*128 + 33], arr[o4*128 + 37], arr_index[o4*128 + 33], arr_index[o4*128 + 37], o4%2==0);
        compareSwap(arr[o4*128 + 34], arr[o4*128 + 38], arr_index[o4*128 + 34], arr_index[o4*128 + 38], o4%2==0);
        compareSwap(arr[o4*128 + 35], arr[o4*128 + 39], arr_index[o4*128 + 35], arr_index[o4*128 + 39], o4%2==0);

        compareSwap(arr[o4*128 + 40], arr[o4*128 + 44], arr_index[o4*128 + 40], arr_index[o4*128 + 44], o4%2==0);
        compareSwap(arr[o4*128 + 41], arr[o4*128 + 45], arr_index[o4*128 + 41], arr_index[o4*128 + 45], o4%2==0);
        compareSwap(arr[o4*128 + 42], arr[o4*128 + 46], arr_index[o4*128 + 42], arr_index[o4*128 + 46], o4%2==0);
        compareSwap(arr[o4*128 + 43], arr[o4*128 + 47], arr_index[o4*128 + 43], arr_index[o4*128 + 47], o4%2==0);

        compareSwap(arr[o4*128 + 48], arr[o4*128 + 52], arr_index[o4*128 + 48], arr_index[o4*128 + 52], o4%2==0);
        compareSwap(arr[o4*128 + 49], arr[o4*128 + 53], arr_index[o4*128 + 49], arr_index[o4*128 + 53], o4%2==0);
        compareSwap(arr[o4*128 + 50], arr[o4*128 + 54], arr_index[o4*128 + 50], arr_index[o4*128 + 54], o4%2==0);
        compareSwap(arr[o4*128 + 51], arr[o4*128 + 55], arr_index[o4*128 + 51], arr_index[o4*128 + 55], o4%2==0);

        compareSwap(arr[o4*128 + 56], arr[o4*128 + 60], arr_index[o4*128 + 56], arr_index[o4*128 + 60], o4%2==0);
        compareSwap(arr[o4*128 + 57], arr[o4*128 + 61], arr_index[o4*128 + 57], arr_index[o4*128 + 61], o4%2==0);
        compareSwap(arr[o4*128 + 58], arr[o4*128 + 62], arr_index[o4*128 + 58], arr_index[o4*128 + 62], o4%2==0);
        compareSwap(arr[o4*128 + 59], arr[o4*128 + 63], arr_index[o4*128 + 59], arr_index[o4*128 + 63], o4%2==0);

        compareSwap(arr[o4*128 + 64], arr[o4*128 + 68], arr_index[o4*128 + 64], arr_index[o4*128 + 68], o4%2==0);
        compareSwap(arr[o4*128 + 65], arr[o4*128 + 69], arr_index[o4*128 + 65], arr_index[o4*128 + 69], o4%2==0);
        compareSwap(arr[o4*128 + 66], arr[o4*128 + 70], arr_index[o4*128 + 66], arr_index[o4*128 + 70], o4%2==0);
        compareSwap(arr[o4*128 + 67], arr[o4*128 + 71], arr_index[o4*128 + 67], arr_index[o4*128 + 71], o4%2==0);

        compareSwap(arr[o4*128 + 72], arr[o4*128 + 76], arr_index[o4*128 + 72], arr_index[o4*128 + 76], o4%2==0);
        compareSwap(arr[o4*128 + 73], arr[o4*128 + 77], arr_index[o4*128 + 73], arr_index[o4*128 + 77], o4%2==0);
        compareSwap(arr[o4*128 + 74], arr[o4*128 + 78], arr_index[o4*128 + 74], arr_index[o4*128 + 78], o4%2==0);
        compareSwap(arr[o4*128 + 75], arr[o4*128 + 79], arr_index[o4*128 + 75], arr_index[o4*128 + 79], o4%2==0);

        compareSwap(arr[o4*128 + 80], arr[o4*128 + 84], arr_index[o4*128 + 80], arr_index[o4*128 + 84], o4%2==0);
        compareSwap(arr[o4*128 + 81], arr[o4*128 + 85], arr_index[o4*128 + 81], arr_index[o4*128 + 85], o4%2==0);
        compareSwap(arr[o4*128 + 82], arr[o4*128 + 86], arr_index[o4*128 + 82], arr_index[o4*128 + 86], o4%2==0);
        compareSwap(arr[o4*128 + 83], arr[o4*128 + 87], arr_index[o4*128 + 83], arr_index[o4*128 + 87], o4%2==0);

        compareSwap(arr[o4*128 + 88], arr[o4*128 + 92], arr_index[o4*128 + 88], arr_index[o4*128 + 92], o4%2==0);
        compareSwap(arr[o4*128 + 89], arr[o4*128 + 93], arr_index[o4*128 + 89], arr_index[o4*128 + 93], o4%2==0);
        compareSwap(arr[o4*128 + 90], arr[o4*128 + 94], arr_index[o4*128 + 90], arr_index[o4*128 + 94], o4%2==0);
        compareSwap(arr[o4*128 + 91], arr[o4*128 + 95], arr_index[o4*128 + 91], arr_index[o4*128 + 95], o4%2==0);

        compareSwap(arr[o4*128 + 96], arr[o4*128 + 100], arr_index[o4*128 + 96], arr_index[o4*128 + 100], o4%2==0);
        compareSwap(arr[o4*128 + 97], arr[o4*128 + 101], arr_index[o4*128 + 97], arr_index[o4*128 + 101], o4%2==0);
        compareSwap(arr[o4*128 + 98], arr[o4*128 + 102], arr_index[o4*128 + 98], arr_index[o4*128 + 102], o4%2==0);
        compareSwap(arr[o4*128 + 99], arr[o4*128 + 103], arr_index[o4*128 + 99], arr_index[o4*128 + 103], o4%2==0);

        compareSwap(arr[o4*128 + 104], arr[o4*128 + 108], arr_index[o4*128 + 104], arr_index[o4*128 + 108], o4%2==0);
        compareSwap(arr[o4*128 + 105], arr[o4*128 + 109], arr_index[o4*128 + 105], arr_index[o4*128 + 109], o4%2==0);
        compareSwap(arr[o4*128 + 106], arr[o4*128 + 110], arr_index[o4*128 + 106], arr_index[o4*128 + 110], o4%2==0);
        compareSwap(arr[o4*128 + 107], arr[o4*128 + 111], arr_index[o4*128 + 107], arr_index[o4*128 + 111], o4%2==0);

        compareSwap(arr[o4*128 + 112], arr[o4*128 + 116], arr_index[o4*128 + 112], arr_index[o4*128 + 116], o4%2==0);
        compareSwap(arr[o4*128 + 113], arr[o4*128 + 117], arr_index[o4*128 + 113], arr_index[o4*128 + 117], o4%2==0);
        compareSwap(arr[o4*128 + 114], arr[o4*128 + 118], arr_index[o4*128 + 114], arr_index[o4*128 + 118], o4%2==0);
        compareSwap(arr[o4*128 + 115], arr[o4*128 + 119], arr_index[o4*128 + 115], arr_index[o4*128 + 119], o4%2==0);

        compareSwap(arr[o4*128 + 120], arr[o4*128 + 124], arr_index[o4*128 + 120], arr_index[o4*128 + 124], o4%2==0);
        compareSwap(arr[o4*128 + 121], arr[o4*128 + 125], arr_index[o4*128 + 121], arr_index[o4*128 + 125], o4%2==0);
        compareSwap(arr[o4*128 + 122], arr[o4*128 + 126], arr_index[o4*128 + 122], arr_index[o4*128 + 126], o4%2==0);
        compareSwap(arr[o4*128 + 123], arr[o4*128 + 127], arr_index[o4*128 + 123], arr_index[o4*128 + 127], o4%2==0);
    }
    for(int o5 = 0; o5 < 1; o5++){

        compareSwap(arr[o5*128], arr[o5*128 + 2], arr_index[o5*128], arr_index[o5*128 + 2], o5%2==0);
        compareSwap(arr[o5*128 + 1], arr[o5*128 + 3], arr_index[o5*128 + 1], arr_index[o5*128 + 3], o5%2==0);
        compareSwap(arr[o5*128 + 4], arr[o5*128 + 6], arr_index[o5*128 + 4], arr_index[o5*128 + 6], o5%2==0);
        compareSwap(arr[o5*128 + 5], arr[o5*128 + 7], arr_index[o5*128 + 5], arr_index[o5*128 + 7], o5%2==0);
        compareSwap(arr[o5*128 + 8], arr[o5*128 + 10], arr_index[o5*128 + 8], arr_index[o5*128 + 10], o5%2==0);
        compareSwap(arr[o5*128 + 9], arr[o5*128 + 11], arr_index[o5*128 + 9], arr_index[o5*128 + 11], o5%2==0);
        compareSwap(arr[o5*128 + 12], arr[o5*128 + 14], arr_index[o5*128 + 12], arr_index[o5*128 + 14], o5%2==0);
        compareSwap(arr[o5*128 + 13], arr[o5*128 + 15], arr_index[o5*128 + 13], arr_index[o5*128 + 15], o5%2==0);
        compareSwap(arr[o5*128 + 16], arr[o5*128 + 18], arr_index[o5*128 + 16], arr_index[o5*128 + 18], o5%2==0);
        compareSwap(arr[o5*128 + 17], arr[o5*128 + 19], arr_index[o5*128 + 17], arr_index[o5*128 + 19], o5%2==0);
        compareSwap(arr[o5*128 + 20], arr[o5*128 + 22], arr_index[o5*128 + 20], arr_index[o5*128 + 22], o5%2==0);
        compareSwap(arr[o5*128 + 21], arr[o5*128 + 23], arr_index[o5*128 + 21], arr_index[o5*128 + 23], o5%2==0);
        compareSwap(arr[o5*128 + 24], arr[o5*128 + 26], arr_index[o5*128 + 24], arr_index[o5*128 + 26], o5%2==0);
        compareSwap(arr[o5*128 + 25], arr[o5*128 + 27], arr_index[o5*128 + 25], arr_index[o5*128 + 27], o5%2==0);
        compareSwap(arr[o5*128 + 28], arr[o5*128 + 30], arr_index[o5*128 + 28], arr_index[o5*128 + 30], o5%2==0);
        compareSwap(arr[o5*128 + 29], arr[o5*128 + 31], arr_index[o5*128 + 29], arr_index[o5*128 + 31], o5%2==0);
        compareSwap(arr[o5*128 + 32], arr[o5*128 + 34], arr_index[o5*128 + 32], arr_index[o5*128 + 34], o5%2==0);
        compareSwap(arr[o5*128 + 33], arr[o5*128 + 35], arr_index[o5*128 + 33], arr_index[o5*128 + 35], o5%2==0);
        compareSwap(arr[o5*128 + 36], arr[o5*128 + 38], arr_index[o5*128 + 36], arr_index[o5*128 + 38], o5%2==0);
        compareSwap(arr[o5*128 + 37], arr[o5*128 + 39], arr_index[o5*128 + 37], arr_index[o5*128 + 39], o5%2==0);
        compareSwap(arr[o5*128 + 40], arr[o5*128 + 42], arr_index[o5*128 + 40], arr_index[o5*128 + 42], o5%2==0);
        compareSwap(arr[o5*128 + 41], arr[o5*128 + 43], arr_index[o5*128 + 41], arr_index[o5*128 + 43], o5%2==0);
        compareSwap(arr[o5*128 + 44], arr[o5*128 + 46], arr_index[o5*128 + 44], arr_index[o5*128 + 46], o5%2==0);
        compareSwap(arr[o5*128 + 45], arr[o5*128 + 47], arr_index[o5*128 + 45], arr_index[o5*128 + 47], o5%2==0);
        compareSwap(arr[o5*128 + 48], arr[o5*128 + 50], arr_index[o5*128 + 48], arr_index[o5*128 + 50], o5%2==0);
        compareSwap(arr[o5*128 + 49], arr[o5*128 + 51], arr_index[o5*128 + 49], arr_index[o5*128 + 51], o5%2==0);
        compareSwap(arr[o5*128 + 52], arr[o5*128 + 54], arr_index[o5*128 + 52], arr_index[o5*128 + 54], o5%2==0);
        compareSwap(arr[o5*128 + 53], arr[o5*128 + 55], arr_index[o5*128 + 53], arr_index[o5*128 + 55], o5%2==0);
        compareSwap(arr[o5*128 + 56], arr[o5*128 + 58], arr_index[o5*128 + 56], arr_index[o5*128 + 58], o5%2==0);
        compareSwap(arr[o5*128 + 57], arr[o5*128 + 59], arr_index[o5*128 + 57], arr_index[o5*128 + 59], o5%2==0);
        compareSwap(arr[o5*128 + 60], arr[o5*128 + 62], arr_index[o5*128 + 60], arr_index[o5*128 + 62], o5%2==0);
        compareSwap(arr[o5*128 + 61], arr[o5*128 + 63], arr_index[o5*128 + 61], arr_index[o5*128 + 63], o5%2==0);
        compareSwap(arr[o5*128 + 64], arr[o5*128 + 66], arr_index[o5*128 + 64], arr_index[o5*128 + 66], o5%2==0);
        compareSwap(arr[o5*128 + 65], arr[o5*128 + 67], arr_index[o5*128 + 65], arr_index[o5*128 + 67], o5%2==0);
        compareSwap(arr[o5*128 + 68], arr[o5*128 + 70], arr_index[o5*128 + 68], arr_index[o5*128 + 70], o5%2==0);
        compareSwap(arr[o5*128 + 69], arr[o5*128 + 71], arr_index[o5*128 + 69], arr_index[o5*128 + 71], o5%2==0);
        compareSwap(arr[o5*128 + 72], arr[o5*128 + 74], arr_index[o5*128 + 72], arr_index[o5*128 + 74], o5%2==0);
        compareSwap(arr[o5*128 + 73], arr[o5*128 + 75], arr_index[o5*128 + 73], arr_index[o5*128 + 75], o5%2==0);
        compareSwap(arr[o5*128 + 76], arr[o5*128 + 78], arr_index[o5*128 + 76], arr_index[o5*128 + 78], o5%2==0);
        compareSwap(arr[o5*128 + 77], arr[o5*128 + 79], arr_index[o5*128 + 77], arr_index[o5*128 + 79], o5%2==0);
        compareSwap(arr[o5*128 + 80], arr[o5*128 + 82], arr_index[o5*128 + 80], arr_index[o5*128 + 82], o5%2==0);
        compareSwap(arr[o5*128 + 81], arr[o5*128 + 83], arr_index[o5*128 + 81], arr_index[o5*128 + 83], o5%2==0);
        compareSwap(arr[o5*128 + 84], arr[o5*128 + 86], arr_index[o5*128 + 84], arr_index[o5*128 + 86], o5%2==0);
        compareSwap(arr[o5*128 + 85], arr[o5*128 + 87], arr_index[o5*128 + 85], arr_index[o5*128 + 87], o5%2==0);
        compareSwap(arr[o5*128 + 88], arr[o5*128 + 90], arr_index[o5*128 + 88], arr_index[o5*128 + 90], o5%2==0);
        compareSwap(arr[o5*128 + 89], arr[o5*128 + 91], arr_index[o5*128 + 89], arr_index[o5*128 + 91], o5%2==0);
        compareSwap(arr[o5*128 + 92], arr[o5*128 + 94], arr_index[o5*128 + 92], arr_index[o5*128 + 94], o5%2==0);
        compareSwap(arr[o5*128 + 93], arr[o5*128 + 95], arr_index[o5*128 + 93], arr_index[o5*128 + 95], o5%2==0);
        compareSwap(arr[o5*128 + 96], arr[o5*128 + 98], arr_index[o5*128 + 96], arr_index[o5*128 + 98], o5%2==0);
        compareSwap(arr[o5*128 + 97], arr[o5*128 + 99], arr_index[o5*128 + 97], arr_index[o5*128 + 99], o5%2==0);
        compareSwap(arr[o5*128 + 100], arr[o5*128 + 102], arr_index[o5*128 + 100], arr_index[o5*128 + 102], o5%2==0);
        compareSwap(arr[o5*128 + 101], arr[o5*128 + 103], arr_index[o5*128 + 101], arr_index[o5*128 + 103], o5%2==0);
        compareSwap(arr[o5*128 + 104], arr[o5*128 + 106], arr_index[o5*128 + 104], arr_index[o5*128 + 106], o5%2==0);
        compareSwap(arr[o5*128 + 105], arr[o5*128 + 107], arr_index[o5*128 + 105], arr_index[o5*128 + 107], o5%2==0);
        compareSwap(arr[o5*128 + 108], arr[o5*128 + 110], arr_index[o5*128 + 108], arr_index[o5*128 + 110], o5%2==0);
        compareSwap(arr[o5*128 + 109], arr[o5*128 + 111], arr_index[o5*128 + 109], arr_index[o5*128 + 111], o5%2==0);
        compareSwap(arr[o5*128 + 112], arr[o5*128 + 114], arr_index[o5*128 + 112], arr_index[o5*128 + 114], o5%2==0);
        compareSwap(arr[o5*128 + 113], arr[o5*128 + 115], arr_index[o5*128 + 113], arr_index[o5*128 + 115], o5%2==0);
        compareSwap(arr[o5*128 + 116], arr[o5*128 + 118], arr_index[o5*128 + 116], arr_index[o5*128 + 118], o5%2==0);
        compareSwap(arr[o5*128 + 117], arr[o5*128 + 119], arr_index[o5*128 + 117], arr_index[o5*128 + 119], o5%2==0);
        compareSwap(arr[o5*128 + 120], arr[o5*128 + 122], arr_index[o5*128 + 120], arr_index[o5*128 + 122], o5%2==0);
        compareSwap(arr[o5*128 + 121], arr[o5*128 + 123], arr_index[o5*128 + 121], arr_index[o5*128 + 123], o5%2==0);
        compareSwap(arr[o5*128 + 124], arr[o5*128 + 126], arr_index[o5*128 + 124], arr_index[o5*128 + 126], o5%2==0);
        compareSwap(arr[o5*128 + 125], arr[o5*128 + 127], arr_index[o5*128 + 125], arr_index[o5*128 + 127], o5%2==0);
    }
    for(int o6 = 0; o6 < 1; o6++){

        compareSwap(arr[o6*128], arr[o6*128 + 1], arr_index[o6*128], arr_index[o6*128 + 1], o6%2==0);
        compareSwap(arr[o6*128 + 2], arr[o6*128 + 3], arr_index[o6*128 + 2], arr_index[o6*128 + 3], o6%2==0);
        compareSwap(arr[o6*128 + 4], arr[o6*128 + 5], arr_index[o6*128 + 4], arr_index[o6*128 + 5], o6%2==0);
        compareSwap(arr[o6*128 + 6], arr[o6*128 + 7], arr_index[o6*128 + 6], arr_index[o6*128 + 7], o6%2==0);
        compareSwap(arr[o6*128 + 8], arr[o6*128 + 9], arr_index[o6*128 + 8], arr_index[o6*128 + 9], o6%2==0);
        compareSwap(arr[o6*128 + 10], arr[o6*128 + 11], arr_index[o6*128 + 10], arr_index[o6*128 + 11], o6%2==0);
        compareSwap(arr[o6*128 + 12], arr[o6*128 + 13], arr_index[o6*128 + 12], arr_index[o6*128 + 13], o6%2==0);
        compareSwap(arr[o6*128 + 14], arr[o6*128 + 15], arr_index[o6*128 + 14], arr_index[o6*128 + 15], o6%2==0);
        compareSwap(arr[o6*128 + 16], arr[o6*128 + 17], arr_index[o6*128 + 16], arr_index[o6*128 + 17], o6%2==0);
        compareSwap(arr[o6*128 + 18], arr[o6*128 + 19], arr_index[o6*128 + 18], arr_index[o6*128 + 19], o6%2==0);
        compareSwap(arr[o6*128 + 20], arr[o6*128 + 21], arr_index[o6*128 + 20], arr_index[o6*128 + 21], o6%2==0);
        compareSwap(arr[o6*128 + 22], arr[o6*128 + 23], arr_index[o6*128 + 22], arr_index[o6*128 + 23], o6%2==0);
        compareSwap(arr[o6*128 + 24], arr[o6*128 + 25], arr_index[o6*128 + 24], arr_index[o6*128 + 25], o6%2==0);
        compareSwap(arr[o6*128 + 26], arr[o6*128 + 27], arr_index[o6*128 + 26], arr_index[o6*128 + 27], o6%2==0);
        compareSwap(arr[o6*128 + 28], arr[o6*128 + 29], arr_index[o6*128 + 28], arr_index[o6*128 + 29], o6%2==0);
        compareSwap(arr[o6*128 + 30], arr[o6*128 + 31], arr_index[o6*128 + 30], arr_index[o6*128 + 31], o6%2==0);
        compareSwap(arr[o6*128 + 32], arr[o6*128 + 33], arr_index[o6*128 + 32], arr_index[o6*128 + 33], o6%2==0);
        compareSwap(arr[o6*128 + 34], arr[o6*128 + 35], arr_index[o6*128 + 34], arr_index[o6*128 + 35], o6%2==0);
        compareSwap(arr[o6*128 + 36], arr[o6*128 + 37], arr_index[o6*128 + 36], arr_index[o6*128 + 37], o6%2==0);
        compareSwap(arr[o6*128 + 38], arr[o6*128 + 39], arr_index[o6*128 + 38], arr_index[o6*128 + 39], o6%2==0);
        compareSwap(arr[o6*128 + 40], arr[o6*128 + 41], arr_index[o6*128 + 40], arr_index[o6*128 + 41], o6%2==0);
        compareSwap(arr[o6*128 + 42], arr[o6*128 + 43], arr_index[o6*128 + 42], arr_index[o6*128 + 43], o6%2==0);
        compareSwap(arr[o6*128 + 44], arr[o6*128 + 45], arr_index[o6*128 + 44], arr_index[o6*128 + 45], o6%2==0);
        compareSwap(arr[o6*128 + 46], arr[o6*128 + 47], arr_index[o6*128 + 46], arr_index[o6*128 + 47], o6%2==0);
        compareSwap(arr[o6*128 + 48], arr[o6*128 + 49], arr_index[o6*128 + 48], arr_index[o6*128 + 49], o6%2==0);
        compareSwap(arr[o6*128 + 50], arr[o6*128 + 51], arr_index[o6*128 + 50], arr_index[o6*128 + 51], o6%2==0);
        compareSwap(arr[o6*128 + 52], arr[o6*128 + 53], arr_index[o6*128 + 52], arr_index[o6*128 + 53], o6%2==0);
        compareSwap(arr[o6*128 + 54], arr[o6*128 + 55], arr_index[o6*128 + 54], arr_index[o6*128 + 55], o6%2==0);
        compareSwap(arr[o6*128 + 56], arr[o6*128 + 57], arr_index[o6*128 + 56], arr_index[o6*128 + 57], o6%2==0);
        compareSwap(arr[o6*128 + 58], arr[o6*128 + 59], arr_index[o6*128 + 58], arr_index[o6*128 + 59], o6%2==0);
        compareSwap(arr[o6*128 + 60], arr[o6*128 + 61], arr_index[o6*128 + 60], arr_index[o6*128 + 61], o6%2==0);
        compareSwap(arr[o6*128 + 62], arr[o6*128 + 63], arr_index[o6*128 + 62], arr_index[o6*128 + 63], o6%2==0);
        compareSwap(arr[o6*128 + 64], arr[o6*128 + 65], arr_index[o6*128 + 64], arr_index[o6*128 + 65], o6%2==0);
        compareSwap(arr[o6*128 + 66], arr[o6*128 + 67], arr_index[o6*128 + 66], arr_index[o6*128 + 67], o6%2==0);
        compareSwap(arr[o6*128 + 68], arr[o6*128 + 69], arr_index[o6*128 + 68], arr_index[o6*128 + 69], o6%2==0);
        compareSwap(arr[o6*128 + 70], arr[o6*128 + 71], arr_index[o6*128 + 70], arr_index[o6*128 + 71], o6%2==0);
        compareSwap(arr[o6*128 + 72], arr[o6*128 + 73], arr_index[o6*128 + 72], arr_index[o6*128 + 73], o6%2==0);
        compareSwap(arr[o6*128 + 74], arr[o6*128 + 75], arr_index[o6*128 + 74], arr_index[o6*128 + 75], o6%2==0);
        compareSwap(arr[o6*128 + 76], arr[o6*128 + 77], arr_index[o6*128 + 76], arr_index[o6*128 + 77], o6%2==0);
        compareSwap(arr[o6*128 + 78], arr[o6*128 + 79], arr_index[o6*128 + 78], arr_index[o6*128 + 79], o6%2==0);
        compareSwap(arr[o6*128 + 80], arr[o6*128 + 81], arr_index[o6*128 + 80], arr_index[o6*128 + 81], o6%2==0);
        compareSwap(arr[o6*128 + 82], arr[o6*128 + 83], arr_index[o6*128 + 82], arr_index[o6*128 + 83], o6%2==0);
        compareSwap(arr[o6*128 + 84], arr[o6*128 + 85], arr_index[o6*128 + 84], arr_index[o6*128 + 85], o6%2==0);
        compareSwap(arr[o6*128 + 86], arr[o6*128 + 87], arr_index[o6*128 + 86], arr_index[o6*128 + 87], o6%2==0);
        compareSwap(arr[o6*128 + 88], arr[o6*128 + 89], arr_index[o6*128 + 88], arr_index[o6*128 + 89], o6%2==0);
        compareSwap(arr[o6*128 + 90], arr[o6*128 + 91], arr_index[o6*128 + 90], arr_index[o6*128 + 91], o6%2==0);
        compareSwap(arr[o6*128 + 92], arr[o6*128 + 93], arr_index[o6*128 + 92], arr_index[o6*128 + 93], o6%2==0);
        compareSwap(arr[o6*128 + 94], arr[o6*128 + 95], arr_index[o6*128 + 94], arr_index[o6*128 + 95], o6%2==0);
        compareSwap(arr[o6*128 + 96], arr[o6*128 + 97], arr_index[o6*128 + 96], arr_index[o6*128 + 97], o6%2==0);
        compareSwap(arr[o6*128 + 98], arr[o6*128 + 99], arr_index[o6*128 + 98], arr_index[o6*128 + 99], o6%2==0);
        compareSwap(arr[o6*128 + 100], arr[o6*128 + 101], arr_index[o6*128 + 100], arr_index[o6*128 + 101], o6%2==0);
        compareSwap(arr[o6*128 + 102], arr[o6*128 + 103], arr_index[o6*128 + 102], arr_index[o6*128 + 103], o6%2==0);
        compareSwap(arr[o6*128 + 104], arr[o6*128 + 105], arr_index[o6*128 + 104], arr_index[o6*128 + 105], o6%2==0);
        compareSwap(arr[o6*128 + 106], arr[o6*128 + 107], arr_index[o6*128 + 106], arr_index[o6*128 + 107], o6%2==0);
        compareSwap(arr[o6*128 + 108], arr[o6*128 + 109], arr_index[o6*128 + 108], arr_index[o6*128 + 109], o6%2==0);
        compareSwap(arr[o6*128 + 110], arr[o6*128 + 111], arr_index[o6*128 + 110], arr_index[o6*128 + 111], o6%2==0);
        compareSwap(arr[o6*128 + 112], arr[o6*128 + 113], arr_index[o6*128 + 112], arr_index[o6*128 + 113], o6%2==0);
        compareSwap(arr[o6*128 + 114], arr[o6*128 + 115], arr_index[o6*128 + 114], arr_index[o6*128 + 115], o6%2==0);
        compareSwap(arr[o6*128 + 116], arr[o6*128 + 117], arr_index[o6*128 + 116], arr_index[o6*128 + 117], o6%2==0);
        compareSwap(arr[o6*128 + 118], arr[o6*128 + 119], arr_index[o6*128 + 118], arr_index[o6*128 + 119], o6%2==0);
        compareSwap(arr[o6*128 + 120], arr[o6*128 + 121], arr_index[o6*128 + 120], arr_index[o6*128 + 121], o6%2==0);
        compareSwap(arr[o6*128 + 122], arr[o6*128 + 123], arr_index[o6*128 + 122], arr_index[o6*128 + 123], o6%2==0);
        compareSwap(arr[o6*128 + 124], arr[o6*128 + 125], arr_index[o6*128 + 124], arr_index[o6*128 + 125], o6%2==0);
        compareSwap(arr[o6*128 + 126], arr[o6*128 + 127], arr_index[o6*128 + 126], arr_index[o6*128 + 127], o6%2==0);
    }
}

//

void bitonicSort64_hw(type_dist_hw* arr, voxel_int* arr_index) 
{
// 对数组进行完全分区以实现并行访问
    //Stage 1
    for(int i = 0; i < 32; i++) {
#pragma HLS UNROLL
        compareSwap(arr[i*2], arr[i*2 + 1], arr_index[i*2], arr_index[i*2 + 1], i%2==0);
    }

    //Stage 2
    for(int j = 0; j < 16; j++) {
#pragma HLS UNROLL
        compareSwap(arr[j*4], arr[j*4 + 2], arr_index[j*4], arr_index[j*4 + 2], j%2==0);
        compareSwap(arr[j*4 + 1], arr[j*4 + 3], arr_index[j*4 + 1], arr_index[j*4 + 3], j%2==0);
    }
    for(int jj = 0; jj < 16; jj++) {
#pragma HLS UNROLL
        compareSwap(arr[jj*4], arr[jj*4 + 1], arr_index[jj*4], arr_index[jj*4 + 1], jj%2==0);
        compareSwap(arr[jj*4 + 2], arr[jj*4 + 3], arr_index[jj*4 + 2], arr_index[jj*4 + 3], jj%2==0);
    }

    //Stage 3
    for(int k = 0; k < 8; k++) {
#pragma HLS UNROLL
        compareSwap(arr[k*8], arr[k*8 + 4], arr_index[k*8], arr_index[k*8 + 4], k%2==0);
        compareSwap(arr[k*8 + 1], arr[k*8 + 5], arr_index[k*8 + 1], arr_index[k*8 + 5], k%2==0);
        compareSwap(arr[k*8 + 2], arr[k*8 + 6], arr_index[k*8 + 2], arr_index[k*8 + 6], k%2==0);
        compareSwap(arr[k*8 + 3], arr[k*8 + 7], arr_index[k*8 + 3], arr_index[k*8 + 7], k%2==0);
    }
    for(int kk = 0; kk < 8; kk++) {
#pragma HLS UNROLL
        compareSwap(arr[kk*8], arr[kk*8 + 2], arr_index[kk*8], arr_index[kk*8 + 2], kk%2==0);
        compareSwap(arr[kk*8 + 1], arr[kk*8 + 3], arr_index[kk*8 + 1], arr_index[kk*8 + 3], kk%2==0);
        compareSwap(arr[kk*8 + 4], arr[kk*8 + 6], arr_index[kk*8 + 4], arr_index[kk*8 + 6], kk%2==0);
        compareSwap(arr[kk*8 + 5], arr[kk*8 + 7], arr_index[kk*8 + 5], arr_index[kk*8 + 7], kk%2==0);
    }
    for(int kkk = 0; kkk < 8; kkk++) {
#pragma HLS UNROLL
        compareSwap(arr[kkk*8], arr[kkk*8 + 1], arr_index[kkk*8], arr_index[kkk*8 + 1], kkk%2==0);
        compareSwap(arr[kkk*8 + 2], arr[kkk*8 + 3], arr_index[kkk*8 + 2], arr_index[kkk*8 + 3], kkk%2==0);
        compareSwap(arr[kkk*8 + 4], arr[kkk*8 + 5], arr_index[kkk*8 + 4], arr_index[kkk*8 + 5], kkk%2==0);
        compareSwap(arr[kkk*8 + 6], arr[kkk*8 + 7], arr_index[kkk*8 + 6], arr_index[kkk*8 + 7], kkk%2==0);
    }

    //Stage 4
    for(int l = 0; l < 4; l++){
#pragma HLS UNROLL
        compareSwap(arr[l*16], arr[l*16 + 8], arr_index[l*16], arr_index[l*16 + 8], l%2==0);
        compareSwap(arr[l*16 + 1], arr[l*16 + 9], arr_index[l*16 + 1], arr_index[l*16 + 9], l%2==0);
        compareSwap(arr[l*16 + 2], arr[l*16 + 10], arr_index[l*16 + 2], arr_index[l*16 + 10], l%2==0);
        compareSwap(arr[l*16 + 3], arr[l*16 + 11], arr_index[l*16 + 3], arr_index[l*16 + 11], l%2==0);
        compareSwap(arr[l*16 + 4], arr[l*16 + 12], arr_index[l*16 + 4], arr_index[l*16 + 12], l%2==0);
        compareSwap(arr[l*16 + 5], arr[l*16 + 13], arr_index[l*16 + 5], arr_index[l*16 + 13], l%2==0);
        compareSwap(arr[l*16 + 6], arr[l*16 + 14], arr_index[l*16 + 6], arr_index[l*16 + 14], l%2==0);
        compareSwap(arr[l*16 + 7], arr[l*16 + 15], arr_index[l*16 + 7], arr_index[l*16 + 15], l%2==0);
    }

    for(int ll = 0; ll < 4; ll++){
#pragma HLS UNROLL
        compareSwap(arr[ll*16], arr[ll*16 + 4], arr_index[ll*16], arr_index[ll*16 + 4], ll%2==0);
        compareSwap(arr[ll*16 + 1], arr[ll*16 + 5], arr_index[ll*16 + 1], arr_index[ll*16 + 5], ll%2==0);
        compareSwap(arr[ll*16 + 2], arr[ll*16 + 6], arr_index[ll*16 + 2], arr_index[ll*16 + 6], ll%2==0);
        compareSwap(arr[ll*16 + 3], arr[ll*16 + 7], arr_index[ll*16 + 3], arr_index[ll*16 + 7], ll%2==0);
        compareSwap(arr[ll*16 + 8], arr[ll*16 + 12], arr_index[ll*16 + 8], arr_index[ll*16 + 12], ll%2==0);
        compareSwap(arr[ll*16 + 9], arr[ll*16 + 13], arr_index[ll*16 + 9], arr_index[ll*16 + 13], ll%2==0);
        compareSwap(arr[ll*16 + 10], arr[ll*16 + 14], arr_index[ll*16 + 10], arr_index[ll*16 + 14], ll%2==0);
        compareSwap(arr[ll*16 + 11], arr[ll*16 + 15], arr_index[ll*16 + 11], arr_index[ll*16 + 15], ll%2==0);
    }

    for (int lll = 0; lll < 4; lll++){
#pragma HLS UNROLL
        compareSwap(arr[lll*16], arr[lll*16 + 2], arr_index[lll*16], arr_index[lll*16 + 2], lll%2==0);
        compareSwap(arr[lll*16 + 1], arr[lll*16 + 3], arr_index[lll*16 + 1], arr_index[lll*16 + 3], lll%2==0);
        compareSwap(arr[lll*16 + 4], arr[lll*16 + 6], arr_index[lll*16 + 4], arr_index[lll*16 + 6], lll%2==0);
        compareSwap(arr[lll*16 + 5], arr[lll*16 + 7], arr_index[lll*16 + 5], arr_index[lll*16 + 7], lll%2==0);
        compareSwap(arr[lll*16 + 8], arr[lll*16 + 10], arr_index[lll*16 + 8], arr_index[lll*16 + 10], lll%2==0);
        compareSwap(arr[lll*16 + 9], arr[lll*16 + 11], arr_index[lll*16 + 9], arr_index[lll*16 + 11], lll%2==0);
        compareSwap(arr[lll*16 + 12], arr[lll*16 + 14], arr_index[lll*16 + 12], arr_index[lll*16 + 14], lll%2==0);
        compareSwap(arr[lll*16 + 13], arr[lll*16 + 15], arr_index[lll*16 + 13], arr_index[lll*16 + 15], lll%2==0);
    }

    for (int llll = 0; llll < 4; llll++){
#pragma HLS UNROLL
        compareSwap(arr[llll*16], arr[llll*16 + 1], arr_index[llll*16], arr_index[llll*16 + 1], llll%2==0);
        compareSwap(arr[llll*16 + 2], arr[llll*16 + 3], arr_index[llll*16 + 2], arr_index[llll*16 + 3], llll%2==0);
        compareSwap(arr[llll*16 + 4], arr[llll*16 + 5], arr_index[llll*16 + 4], arr_index[llll*16 + 5], llll%2==0);
        compareSwap(arr[llll*16 + 6], arr[llll*16 + 7], arr_index[llll*16 + 6], arr_index[llll*16 + 7], llll%2==0);
        compareSwap(arr[llll*16 + 8], arr[llll*16 + 9], arr_index[llll*16 + 8], arr_index[llll*16 + 9], llll%2==0);
        compareSwap(arr[llll*16 + 10], arr[llll*16 + 11], arr_index[llll*16 + 10], arr_index[llll*16 + 11], llll%2==0);
        compareSwap(arr[llll*16 + 12], arr[llll*16 + 13], arr_index[llll*16 + 12], arr_index[llll*16 + 13], llll%2==0);
        compareSwap(arr[llll*16 + 14], arr[llll*16 + 15], arr_index[llll*16 + 14], arr_index[llll*16 + 15], llll%2==0);
    }

    //Stage 5

    for(int m = 0; m < 2; m++){
#pragma HLS UNROLL
        compareSwap(arr[m*32], arr[m*32 + 16], arr_index[m*32], arr_index[m*32 + 16], m%2==0);
        compareSwap(arr[m*32 + 1], arr[m*32 + 17], arr_index[m*32 + 1], arr_index[m*32 + 17], m%2==0);
        compareSwap(arr[m*32 + 2], arr[m*32 + 18], arr_index[m*32 + 2], arr_index[m*32 + 18], m%2==0);
        compareSwap(arr[m*32 + 3], arr[m*32 + 19], arr_index[m*32 + 3], arr_index[m*32 + 19], m%2==0);
        compareSwap(arr[m*32 + 4], arr[m*32 + 20], arr_index[m*32 + 4], arr_index[m*32 + 20], m%2==0);
        compareSwap(arr[m*32 + 5], arr[m*32 + 21], arr_index[m*32 + 5], arr_index[m*32 + 21], m%2==0);
        compareSwap(arr[m*32 + 6], arr[m*32 + 22], arr_index[m*32 + 6], arr_index[m*32 + 22], m%2==0);
        compareSwap(arr[m*32 + 7], arr[m*32 + 23], arr_index[m*32 + 7], arr_index[m*32 + 23], m%2==0);
        compareSwap(arr[m*32 + 8], arr[m*32 + 24], arr_index[m*32 + 8], arr_index[m*32 + 24], m%2==0);
        compareSwap(arr[m*32 + 9], arr[m*32 + 25], arr_index[m*32 + 9], arr_index[m*32 + 25], m%2==0);
        compareSwap(arr[m*32 + 10], arr[m*32 + 26], arr_index[m*32 + 10], arr_index[m*32 + 26], m%2==0);
        compareSwap(arr[m*32 + 11], arr[m*32 + 27], arr_index[m*32 + 11], arr_index[m*32 + 27], m%2==0);
        compareSwap(arr[m*32 + 12], arr[m*32 + 28], arr_index[m*32 + 12], arr_index[m*32 + 28], m%2==0);
        compareSwap(arr[m*32 + 13], arr[m*32 + 29], arr_index[m*32 + 13], arr_index[m*32 + 29], m%2==0);
        compareSwap(arr[m*32 + 14], arr[m*32 + 30], arr_index[m*32 + 14], arr_index[m*32 + 30], m%2==0);
        compareSwap(arr[m*32 + 15], arr[m*32 + 31], arr_index[m*32 + 15], arr_index[m*32 + 31], m%2==0);
    }

    for(int mm = 0; mm < 2; mm++){
#pragma HLS UNROLL
        compareSwap(arr[mm*32], arr[mm*32 + 8], arr_index[mm*32], arr_index[mm*32 + 8], mm%2==0);
        compareSwap(arr[mm*32 + 1], arr[mm*32 + 9], arr_index[mm*32 + 1], arr_index[mm*32 + 9], mm%2==0);
        compareSwap(arr[mm*32 + 2], arr[mm*32 + 10], arr_index[mm*32 + 2], arr_index[mm*32 + 10], mm%2==0);
        compareSwap(arr[mm*32 + 3], arr[mm*32 + 11], arr_index[mm*32 + 3], arr_index[mm*32 + 11], mm%2==0);
        compareSwap(arr[mm*32 + 4], arr[mm*32 + 12], arr_index[mm*32 + 4], arr_index[mm*32 + 12], mm%2==0);
        compareSwap(arr[mm*32 + 5], arr[mm*32 + 13], arr_index[mm*32 + 5], arr_index[mm*32 + 13], mm%2==0);
        compareSwap(arr[mm*32 + 6], arr[mm*32 + 14], arr_index[mm*32 + 6], arr_index[mm*32 + 14], mm%2==0);
        compareSwap(arr[mm*32 + 7], arr[mm*32 + 15], arr_index[mm*32 + 7], arr_index[mm*32 + 15], mm%2==0);
        compareSwap(arr[mm*32 + 16], arr[mm*32 + 24], arr_index[mm*32 + 16], arr_index[mm*32 + 24], mm%2==0);
        compareSwap(arr[mm*32 + 17], arr[mm*32 + 25], arr_index[mm*32 + 17], arr_index[mm*32 + 25], mm%2==0);
        compareSwap(arr[mm*32 + 18], arr[mm*32 + 26], arr_index[mm*32 + 18], arr_index[mm*32 + 26], mm%2==0);
        compareSwap(arr[mm*32 + 19], arr[mm*32 + 27], arr_index[mm*32 + 19], arr_index[mm*32 + 27], mm%2==0);
        compareSwap(arr[mm*32 + 20], arr[mm*32 + 28], arr_index[mm*32 + 20], arr_index[mm*32 + 28], mm%2==0);
        compareSwap(arr[mm*32 + 21], arr[mm*32 + 29], arr_index[mm*32 + 21], arr_index[mm*32 + 29], mm%2==0);
        compareSwap(arr[mm*32 + 22], arr[mm*32 + 30], arr_index[mm*32 + 22], arr_index[mm*32 + 30], mm%2==0);
        compareSwap(arr[mm*32 + 23], arr[mm*32 + 31], arr_index[mm*32 + 23], arr_index[mm*32 + 31], mm%2==0);
    }

    for(int mmm = 0; mmm < 2; mmm++){
#pragma HLS UNROLL
        compareSwap(arr[mmm*32], arr[mmm*32 + 4], arr_index[mmm*32], arr_index[mmm*32 + 4], mmm%2==0);
        compareSwap(arr[mmm*32 + 1], arr[mmm*32 + 5], arr_index[mmm*32 + 1], arr_index[mmm*32 + 5], mmm%2==0);
        compareSwap(arr[mmm*32 + 2], arr[mmm*32 + 6], arr_index[mmm*32 + 2], arr_index[mmm*32 + 6], mmm%2==0);
        compareSwap(arr[mmm*32 + 3], arr[mmm*32 + 7], arr_index[mmm*32 + 3], arr_index[mmm*32 + 7], mmm%2==0);
        compareSwap(arr[mmm*32 + 8], arr[mmm*32 + 12], arr_index[mmm*32 + 8], arr_index[mmm*32 + 12], mmm%2==0);
        compareSwap(arr[mmm*32 + 9], arr[mmm*32 + 13], arr_index[mmm*32 + 9], arr_index[mmm*32 + 13], mmm%2==0);
        compareSwap(arr[mmm*32 + 10], arr[mmm*32 + 14], arr_index[mmm*32 + 10], arr_index[mmm*32 + 14], mmm%2==0);
        compareSwap(arr[mmm*32 + 11], arr[mmm*32 + 15], arr_index[mmm*32 + 11], arr_index[mmm*32 + 15], mmm%2==0);
        compareSwap(arr[mmm*32 + 16], arr[mmm*32 + 20], arr_index[mmm*32 + 16], arr_index[mmm*32 + 20], mmm%2==0);
        compareSwap(arr[mmm*32 + 17], arr[mmm*32 + 21], arr_index[mmm*32 + 17], arr_index[mmm*32 + 21], mmm%2==0);
        compareSwap(arr[mmm*32 + 18], arr[mmm*32 + 22], arr_index[mmm*32 + 18], arr_index[mmm*32 + 22], mmm%2==0);
        compareSwap(arr[mmm*32 + 19], arr[mmm*32 + 23], arr_index[mmm*32 + 19], arr_index[mmm*32 + 23], mmm%2==0);
        compareSwap(arr[mmm*32 + 24], arr[mmm*32 + 28], arr_index[mmm*32 + 24], arr_index[mmm*32 + 28], mmm%2==0);
        compareSwap(arr[mmm*32 + 25], arr[mmm*32 + 29], arr_index[mmm*32 + 25], arr_index[mmm*32 + 29], mmm%2==0);
        compareSwap(arr[mmm*32 + 26], arr[mmm*32 + 30], arr_index[mmm*32 + 26], arr_index[mmm*32 + 30], mmm%2==0);
        compareSwap(arr[mmm*32 + 27], arr[mmm*32 + 31], arr_index[mmm*32 + 27], arr_index[mmm*32 + 31], mmm%2==0);
    }

    for(int mmmm = 0; mmmm < 2; mmmm ++){
#pragma HLS UNROLL
        compareSwap(arr[mmmm*32], arr[mmmm*32 + 2], arr_index[mmmm*32], arr_index[mmmm*32 + 2], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 1], arr[mmmm*32 + 3], arr_index[mmmm*32 + 1], arr_index[mmmm*32 + 3], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 4], arr[mmmm*32 + 6], arr_index[mmmm*32 + 4], arr_index[mmmm*32 + 6], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 5], arr[mmmm*32 + 7], arr_index[mmmm*32 + 5], arr_index[mmmm*32 + 7], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 8], arr[mmmm*32 + 10], arr_index[mmmm*32 + 8], arr_index[mmmm*32 + 10], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 9], arr[mmmm*32 + 11], arr_index[mmmm*32 + 9], arr_index[mmmm*32 + 11], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 12], arr[mmmm*32 + 14], arr_index[mmmm*32 + 12], arr_index[mmmm*32 + 14], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 13], arr[mmmm*32 + 15], arr_index[mmmm*32 + 13], arr_index[mmmm*32 + 15], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 16], arr[mmmm*32 + 18], arr_index[mmmm*32 + 16], arr_index[mmmm*32 + 18], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 17], arr[mmmm*32 + 19], arr_index[mmmm*32 + 17], arr_index[mmmm*32 + 19], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 20], arr[mmmm*32 + 22], arr_index[mmmm*32 + 20], arr_index[mmmm*32 + 22], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 21], arr[mmmm*32 + 23], arr_index[mmmm*32 + 21], arr_index[mmmm*32 + 23], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 24], arr[mmmm*32 + 26], arr_index[mmmm*32 + 24], arr_index[mmmm*32 + 26], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 25], arr[mmmm*32 + 27], arr_index[mmmm*32 + 25], arr_index[mmmm*32 + 27], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 28], arr[mmmm*32 + 30], arr_index[mmmm*32 + 28], arr_index[mmmm*32 + 30], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 29], arr[mmmm*32 + 31], arr_index[mmmm*32 + 29], arr_index[mmmm*32 + 31], mmmm%2==0);
    }

    for(int m5 = 0; m5 < 2; m5++){
#pragma HLS UNROLL
        compareSwap(arr[m5*32], arr[m5*32 + 1], arr_index[m5*32], arr_index[m5*32 + 1], m5%2==0);
        compareSwap(arr[m5*32 + 2], arr[m5*32 + 3], arr_index[m5*32 + 2], arr_index[m5*32 + 3], m5%2==0);
        compareSwap(arr[m5*32 + 4], arr[m5*32 + 5], arr_index[m5*32 + 4], arr_index[m5*32 + 5], m5%2==0);
        compareSwap(arr[m5*32 + 6], arr[m5*32 + 7], arr_index[m5*32 + 6], arr_index[m5*32 + 7], m5%2==0);
        compareSwap(arr[m5*32 + 8], arr[m5*32 + 9], arr_index[m5*32 + 8], arr_index[m5*32 + 9], m5%2==0);
        compareSwap(arr[m5*32 + 10], arr[m5*32 + 11], arr_index[m5*32 + 10], arr_index[m5*32 + 11], m5%2==0);
        compareSwap(arr[m5*32 + 12], arr[m5*32 + 13], arr_index[m5*32 + 12], arr_index[m5*32 + 13], m5%2==0);
        compareSwap(arr[m5*32 + 14], arr[m5*32 + 15], arr_index[m5*32 + 14], arr_index[m5*32 + 15], m5%2==0);
        compareSwap(arr[m5*32 + 16], arr[m5*32 + 17], arr_index[m5*32 + 16], arr_index[m5*32 + 17], m5%2==0);
        compareSwap(arr[m5*32 + 18], arr[m5*32 + 19], arr_index[m5*32 + 18], arr_index[m5*32 + 19], m5%2==0);
        compareSwap(arr[m5*32 + 20], arr[m5*32 + 21], arr_index[m5*32 + 20], arr_index[m5*32 + 21], m5%2==0);
        compareSwap(arr[m5*32 + 22], arr[m5*32 + 23], arr_index[m5*32 + 22], arr_index[m5*32 + 23], m5%2==0);
        compareSwap(arr[m5*32 + 24], arr[m5*32 + 25], arr_index[m5*32 + 24], arr_index[m5*32 + 25], m5%2==0);
        compareSwap(arr[m5*32 + 26], arr[m5*32 + 27], arr_index[m5*32 + 26], arr_index[m5*32 + 27], m5%2==0);
        compareSwap(arr[m5*32 + 28], arr[m5*32 + 29], arr_index[m5*32 + 28], arr_index[m5*32 + 29], m5%2==0);
        compareSwap(arr[m5*32 + 30], arr[m5*32 + 31], arr_index[m5*32 + 30], arr_index[m5*32 + 31], m5%2==0);
    }

    //Stage 6

//     for(int n = 0; n < 1; n++){
// #pragma HLS UNROLL
//         compareSwap(arr[n*64], arr[n*64 + 32], arr_index[n*64], arr_index[n*64 + 32], n%2==0);
//         compareSwap(arr[n*64 + 1], arr[n*64 + 33], arr_index[n*64 + 1], arr_index[n*64 + 33], n%2==0);
//         compareSwap(arr[n*64 + 2], arr[n*64 + 34], arr_index[n*64 + 2], arr_index[n*64 + 34], n%2==0);
//         compareSwap(arr[n*64 + 3], arr[n*64 + 35], arr_index[n*64 + 3], arr_index[n*64 + 35], n%2==0);
//         compareSwap(arr[n*64 + 4], arr[n*64 + 36], arr_index[n*64 + 4], arr_index[n*64 + 36], n%2==0);
//         compareSwap(arr[n*64 + 5], arr[n*64 + 37], arr_index[n*64 + 5], arr_index[n*64 + 37], n%2==0);
//         compareSwap(arr[n*64 + 6], arr[n*64 + 38], arr_index[n*64 + 6], arr_index[n*64 + 38], n%2==0);
//         compareSwap(arr[n*64 + 7], arr[n*64 + 39], arr_index[n*64 + 7], arr_index[n*64 + 39], n%2==0);
//         compareSwap(arr[n*64 + 8], arr[n*64 + 40], arr_index[n*64 + 8], arr_index[n*64 + 40], n%2==0);
//         compareSwap(arr[n*64 + 9], arr[n*64 + 41], arr_index[n*64 + 9], arr_index[n*64 + 41], n%2==0);
//         compareSwap(arr[n*64 + 10], arr[n*64 + 42], arr_index[n*64 + 10], arr_index[n*64 + 42], n%2==0);
//         compareSwap(arr[n*64 + 11], arr[n*64 + 43], arr_index[n*64 + 11], arr_index[n*64 + 43], n%2==0);
//         compareSwap(arr[n*64 + 12], arr[n*64 + 44], arr_index[n*64 + 12], arr_index[n*64 + 44], n%2==0);
//         compareSwap(arr[n*64 + 13], arr[n*64 + 45], arr_index[n*64 + 13], arr_index[n*64 + 45], n%2==0);
//         compareSwap(arr[n*64 + 14], arr[n*64 + 46], arr_index[n*64 + 14], arr_index[n*64 + 46], n%2==0);
//         compareSwap(arr[n*64 + 15], arr[n*64 + 47], arr_index[n*64 + 15], arr_index[n*64 + 47], n%2==0);
//         compareSwap(arr[n*64 + 16], arr[n*64 + 48], arr_index[n*64 + 16], arr_index[n*64 + 48], n%2==0);
//         compareSwap(arr[n*64 + 17], arr[n*64 + 49], arr_index[n*64 + 17], arr_index[n*64 + 49], n%2==0);
//         compareSwap(arr[n*64 + 18], arr[n*64 + 50], arr_index[n*64 + 18], arr_index[n*64 + 50], n%2==0);
//         compareSwap(arr[n*64 + 19], arr[n*64 + 51], arr_index[n*64 + 19], arr_index[n*64 + 51], n%2==0);
//         compareSwap(arr[n*64 + 20], arr[n*64 + 52], arr_index[n*64 + 20], arr_index[n*64 + 52], n%2==0);
//         compareSwap(arr[n*64 + 21], arr[n*64 + 53], arr_index[n*64 + 21], arr_index[n*64 + 53], n%2==0);
//         compareSwap(arr[n*64 + 22], arr[n*64 + 54], arr_index[n*64 + 22], arr_index[n*64 + 54], n%2==0);
//         compareSwap(arr[n*64 + 23], arr[n*64 + 55], arr_index[n*64 + 23], arr_index[n*64 + 55], n%2==0);
//         compareSwap(arr[n*64 + 24], arr[n*64 + 56], arr_index[n*64 + 24], arr_index[n*64 + 56], n%2==0);
//         compareSwap(arr[n*64 + 25], arr[n*64 + 57], arr_index[n*64 + 25], arr_index[n*64 + 57], n%2==0);
//         compareSwap(arr[n*64 + 26], arr[n*64 + 58], arr_index[n*64 + 26], arr_index[n*64 + 58], n%2==0);
//         compareSwap(arr[n*64 + 27], arr[n*64 + 59], arr_index[n*64 + 27], arr_index[n*64 + 59], n%2==0);
//         compareSwap(arr[n*64 + 28], arr[n*64 + 60], arr_index[n*64 + 28], arr_index[n*64 + 60], n%2==0);
//         compareSwap(arr[n*64 + 29], arr[n*64 + 61], arr_index[n*64 + 29], arr_index[n*64 + 61], n%2==0);
//         compareSwap(arr[n*64 + 30], arr[n*64 + 62], arr_index[n*64 + 30], arr_index[n*64 + 62], n%2==0);
//         compareSwap(arr[n*64 + 31], arr[n*64 + 63], arr_index[n*64 + 31], arr_index[n*64 + 63], n%2==0);
//     }

//     for(int n2 = 0; n2 < 1; n2++){
// #pragma HLS UNROLL
//         compareSwap(arr[n2*64], arr[n2*64 + 16], arr_index[n2*64], arr_index[n2*64 + 16], n2%2==0);
//         compareSwap(arr[n2*64 + 1], arr[n2*64 + 17], arr_index[n2*64 + 1], arr_index[n2*64 + 17], n2%2==0);
//         compareSwap(arr[n2*64 + 2], arr[n2*64 + 18], arr_index[n2*64 + 2], arr_index[n2*64 + 18], n2%2==0);
//         compareSwap(arr[n2*64 + 3], arr[n2*64 + 19], arr_index[n2*64 + 3], arr_index[n2*64 + 19], n2%2==0);
//         compareSwap(arr[n2*64 + 4], arr[n2*64 + 20], arr_index[n2*64 + 4], arr_index[n2*64 + 20], n2%2==0);
//         compareSwap(arr[n2*64 + 5], arr[n2*64 + 21], arr_index[n2*64 + 5], arr_index[n2*64 + 21], n2%2==0);
//         compareSwap(arr[n2*64 + 6], arr[n2*64 + 22], arr_index[n2*64 + 6], arr_index[n2*64 + 22], n2%2==0);
//         compareSwap(arr[n2*64 + 7], arr[n2*64 + 23], arr_index[n2*64 + 7], arr_index[n2*64 + 23], n2%2==0);
//         compareSwap(arr[n2*64 + 8], arr[n2*64 + 24], arr_index[n2*64 + 8], arr_index[n2*64 + 24], n2%2==0);
//         compareSwap(arr[n2*64 + 9], arr[n2*64 + 25], arr_index[n2*64 + 9], arr_index[n2*64 + 25], n2%2==0);
//         compareSwap(arr[n2*64 + 10], arr[n2*64 + 26], arr_index[n2*64 + 10], arr_index[n2*64 + 26], n2%2==0);
//         compareSwap(arr[n2*64 + 11], arr[n2*64 + 27], arr_index[n2*64 + 11], arr_index[n2*64 + 27], n2%2==0);
//         compareSwap(arr[n2*64 + 12], arr[n2*64 + 28], arr_index[n2*64 + 12], arr_index[n2*64 + 28], n2%2==0);
//         compareSwap(arr[n2*64 + 13], arr[n2*64 + 29], arr_index[n2*64 + 13], arr_index[n2*64 + 29], n2%2==0);
//         compareSwap(arr[n2*64 + 14], arr[n2*64 + 30], arr_index[n2*64 + 14], arr_index[n2*64 + 30], n2%2==0);
//         compareSwap(arr[n2*64 + 15], arr[n2*64 + 31], arr_index[n2*64 + 15], arr_index[n2*64 + 31], n2%2==0);
        
//         compareSwap(arr[n2*64 + 32], arr[n2*64 + 48], arr_index[n2*64 + 32], arr_index[n2*64 + 48], n2%2==0);
//         compareSwap(arr[n2*64 + 33], arr[n2*64 + 49], arr_index[n2*64 + 33], arr_index[n2*64 + 49], n2%2==0);
//         compareSwap(arr[n2*64 + 34], arr[n2*64 + 50], arr_index[n2*64 + 34], arr_index[n2*64 + 50], n2%2==0);
//         compareSwap(arr[n2*64 + 35], arr[n2*64 + 51], arr_index[n2*64 + 35], arr_index[n2*64 + 51], n2%2==0);
//         compareSwap(arr[n2*64 + 36], arr[n2*64 + 52], arr_index[n2*64 + 36], arr_index[n2*64 + 52], n2%2==0);
//         compareSwap(arr[n2*64 + 37], arr[n2*64 + 53], arr_index[n2*64 + 37], arr_index[n2*64 + 53], n2%2==0);
//         compareSwap(arr[n2*64 + 38], arr[n2*64 + 54], arr_index[n2*64 + 38], arr_index[n2*64 + 54], n2%2==0);
//         compareSwap(arr[n2*64 + 39], arr[n2*64 + 55], arr_index[n2*64 + 39], arr_index[n2*64 + 55], n2%2==0);
//         compareSwap(arr[n2*64 + 40], arr[n2*64 + 56], arr_index[n2*64 + 40], arr_index[n2*64 + 56], n2%2==0);
//         compareSwap(arr[n2*64 + 41], arr[n2*64 + 57], arr_index[n2*64 + 41], arr_index[n2*64 + 57], n2%2==0);
//         compareSwap(arr[n2*64 + 42], arr[n2*64 + 58], arr_index[n2*64 + 42], arr_index[n2*64 + 58], n2%2==0);
//         compareSwap(arr[n2*64 + 43], arr[n2*64 + 59], arr_index[n2*64 + 43], arr_index[n2*64 + 59], n2%2==0);
//         compareSwap(arr[n2*64 + 44], arr[n2*64 + 60], arr_index[n2*64 + 44], arr_index[n2*64 + 60], n2%2==0);
//         compareSwap(arr[n2*64 + 45], arr[n2*64 + 61], arr_index[n2*64 + 45], arr_index[n2*64 + 61], n2%2==0);
//         compareSwap(arr[n2*64 + 46], arr[n2*64 + 62], arr_index[n2*64 + 46], arr_index[n2*64 + 62], n2%2==0);
//         compareSwap(arr[n2*64 + 47], arr[n2*64 + 63], arr_index[n2*64 + 47], arr_index[n2*64 + 63], n2%2==0);
//     }

//     for(int n3 = 0; n3 < 1; n3++){
// #pragma HLS UNROLL
//         compareSwap(arr[n3*64], arr[n3*64 + 8], arr_index[n3*64], arr_index[n3*64 + 8], n3%2==0);
//         compareSwap(arr[n3*64 + 1], arr[n3*64 + 9], arr_index[n3*64 + 1], arr_index[n3*64 + 9], n3%2==0);
//         compareSwap(arr[n3*64 + 2], arr[n3*64 + 10], arr_index[n3*64 + 2], arr_index[n3*64 + 10], n3%2==0);
//         compareSwap(arr[n3*64 + 3], arr[n3*64 + 11], arr_index[n3*64 + 3], arr_index[n3*64 + 11], n3%2==0);
//         compareSwap(arr[n3*64 + 4], arr[n3*64 + 12], arr_index[n3*64 + 4], arr_index[n3*64 + 12], n3%2==0);
//         compareSwap(arr[n3*64 + 5], arr[n3*64 + 13], arr_index[n3*64 + 5], arr_index[n3*64 + 13], n3%2==0);
//         compareSwap(arr[n3*64 + 6], arr[n3*64 + 14], arr_index[n3*64 + 6], arr_index[n3*64 + 14], n3%2==0);
//         compareSwap(arr[n3*64 + 7], arr[n3*64 + 15], arr_index[n3*64 + 7], arr_index[n3*64 + 15], n3%2==0);

//         compareSwap(arr[n3*64 + 16], arr[n3*64 + 24], arr_index[n3*64 + 16], arr_index[n3*64 + 24], n3%2==0);
//         compareSwap(arr[n3*64 + 17], arr[n3*64 + 25], arr_index[n3*64 + 17], arr_index[n3*64 + 25], n3%2==0);
//         compareSwap(arr[n3*64 + 18], arr[n3*64 + 26], arr_index[n3*64 + 18], arr_index[n3*64 + 26], n3%2==0);
//         compareSwap(arr[n3*64 + 19], arr[n3*64 + 27], arr_index[n3*64 + 19], arr_index[n3*64 + 27], n3%2==0);
//         compareSwap(arr[n3*64 + 20], arr[n3*64 + 28], arr_index[n3*64 + 20], arr_index[n3*64 + 28], n3%2==0);
//         compareSwap(arr[n3*64 + 21], arr[n3*64 + 29], arr_index[n3*64 + 21], arr_index[n3*64 + 29], n3%2==0);
//         compareSwap(arr[n3*64 + 22], arr[n3*64 + 30], arr_index[n3*64 + 22], arr_index[n3*64 + 30], n3%2==0);
//         compareSwap(arr[n3*64 + 23], arr[n3*64 + 31], arr_index[n3*64 + 23], arr_index[n3*64 + 31], n3%2==0);

//         compareSwap(arr[n3*64 + 32], arr[n3*64 + 40], arr_index[n3*64 + 32], arr_index[n3*64 + 40], n3%2==0);
//         compareSwap(arr[n3*64 + 33], arr[n3*64 + 41], arr_index[n3*64 + 33], arr_index[n3*64 + 41], n3%2==0);
//         compareSwap(arr[n3*64 + 34], arr[n3*64 + 42], arr_index[n3*64 + 34], arr_index[n3*64 + 42], n3%2==0);
//         compareSwap(arr[n3*64 + 35], arr[n3*64 + 43], arr_index[n3*64 + 35], arr_index[n3*64 + 43], n3%2==0);
//         compareSwap(arr[n3*64 + 36], arr[n3*64 + 44], arr_index[n3*64 + 36], arr_index[n3*64 + 44], n3%2==0);
//         compareSwap(arr[n3*64 + 37], arr[n3*64 + 45], arr_index[n3*64 + 37], arr_index[n3*64 + 45], n3%2==0);
//         compareSwap(arr[n3*64 + 38], arr[n3*64 + 46], arr_index[n3*64 + 38], arr_index[n3*64 + 46], n3%2==0);
//         compareSwap(arr[n3*64 + 39], arr[n3*64 + 47], arr_index[n3*64 + 39], arr_index[n3*64 + 47], n3%2==0);
        
//         compareSwap(arr[n3*64 + 48], arr[n3*64 + 56], arr_index[n3*64 + 48], arr_index[n3*64 + 56], n3%2==0);
//         compareSwap(arr[n3*64 + 49], arr[n3*64 + 57], arr_index[n3*64 + 49], arr_index[n3*64 + 57], n3%2==0);
//         compareSwap(arr[n3*64 + 50], arr[n3*64 + 58], arr_index[n3*64 + 50], arr_index[n3*64 + 58], n3%2==0);
//         compareSwap(arr[n3*64 + 51], arr[n3*64 + 59], arr_index[n3*64 + 51], arr_index[n3*64 + 59], n3%2==0);
//         compareSwap(arr[n3*64 + 52], arr[n3*64 + 60], arr_index[n3*64 + 52], arr_index[n3*64 + 60], n3%2==0);
//         compareSwap(arr[n3*64 + 53], arr[n3*64 + 61], arr_index[n3*64 + 53], arr_index[n3*64 + 61], n3%2==0);
//         compareSwap(arr[n3*64 + 54], arr[n3*64 + 62], arr_index[n3*64 + 54], arr_index[n3*64 + 62], n3%2==0);
//         compareSwap(arr[n3*64 + 55], arr[n3*64 + 63], arr_index[n3*64 + 55], arr_index[n3*64 + 63], n3%2==0);
//     }

//     for(int n4 = 0; n4 < 1; n4++){
// #pragma HLS UNROLL
//         compareSwap(arr[n4*64], arr[n4*64 + 4], arr_index[n4*64], arr_index[n4*64 + 4], n4%2==0);
//         compareSwap(arr[n4*64 + 1], arr[n4*64 + 5], arr_index[n4*64 + 1], arr_index[n4*64 + 5], n4%2==0);
//         compareSwap(arr[n4*64 + 2], arr[n4*64 + 6], arr_index[n4*64 + 2], arr_index[n4*64 + 6], n4%2==0);
//         compareSwap(arr[n4*64 + 3], arr[n4*64 + 7], arr_index[n4*64 + 3], arr_index[n4*64 + 7], n4%2==0);
//         compareSwap(arr[n4*64 + 8], arr[n4*64 + 12], arr_index[n4*64 + 8], arr_index[n4*64 + 12], n4%2==0);
//         compareSwap(arr[n4*64 + 9], arr[n4*64 + 13], arr_index[n4*64 + 9], arr_index[n4*64 + 13], n4%2==0);
//         compareSwap(arr[n4*64 + 10], arr[n4*64 + 14], arr_index[n4*64 + 10], arr_index[n4*64 + 14], n4%2==0);
//         compareSwap(arr[n4*64 + 11], arr[n4*64 + 15], arr_index[n4*64 + 11], arr_index[n4*64 + 15], n4%2==0);
//         compareSwap(arr[n4*64 + 16], arr[n4*64 + 20], arr_index[n4*64 + 16], arr_index[n4*64 + 20], n4%2==0);
//         compareSwap(arr[n4*64 + 17], arr[n4*64 + 21], arr_index[n4*64 + 17], arr_index[n4*64 + 21], n4%2==0);
//         compareSwap(arr[n4*64 + 18], arr[n4*64 + 22], arr_index[n4*64 + 18], arr_index[n4*64 + 22], n4%2==0);
//         compareSwap(arr[n4*64 + 19], arr[n4*64 + 23], arr_index[n4*64 + 19], arr_index[n4*64 + 23], n4%2==0);
//         compareSwap(arr[n4*64 + 24], arr[n4*64 + 28], arr_index[n4*64 + 24], arr_index[n4*64 + 28], n4%2==0);
//         compareSwap(arr[n4*64 + 25], arr[n4*64 + 29], arr_index[n4*64 + 25], arr_index[n4*64 + 29], n4%2==0);
//         compareSwap(arr[n4*64 + 26], arr[n4*64 + 30], arr_index[n4*64 + 26], arr_index[n4*64 + 30], n4%2==0);
//         compareSwap(arr[n4*64 + 27], arr[n4*64 + 31], arr_index[n4*64 + 27], arr_index[n4*64 + 31], n4%2==0);
//         compareSwap(arr[n4*64 + 32], arr[n4*64 + 36], arr_index[n4*64 + 32], arr_index[n4*64 + 36], n4%2==0);
//         compareSwap(arr[n4*64 + 33], arr[n4*64 + 37], arr_index[n4*64 + 33], arr_index[n4*64 + 37], n4%2==0);
//         compareSwap(arr[n4*64 + 34], arr[n4*64 + 38], arr_index[n4*64 + 34], arr_index[n4*64 + 38], n4%2==0);
//         compareSwap(arr[n4*64 + 35], arr[n4*64 + 39], arr_index[n4*64 + 35], arr_index[n4*64 + 39], n4%2==0);

//         compareSwap(arr[n4*64 + 40], arr[n4*64 + 44], arr_index[n4*64 + 40], arr_index[n4*64 + 44], n4%2==0);
//         compareSwap(arr[n4*64 + 41], arr[n4*64 + 45], arr_index[n4*64 + 41], arr_index[n4*64 + 45], n4%2==0);
//         compareSwap(arr[n4*64 + 42], arr[n4*64 + 46], arr_index[n4*64 + 42], arr_index[n4*64 + 46], n4%2==0);
//         compareSwap(arr[n4*64 + 43], arr[n4*64 + 47], arr_index[n4*64 + 43], arr_index[n4*64 + 47], n4%2==0);

//         compareSwap(arr[n4*64 + 48], arr[n4*64 + 52], arr_index[n4*64 + 48], arr_index[n4*64 + 52], n4%2==0);
//         compareSwap(arr[n4*64 + 49], arr[n4*64 + 53], arr_index[n4*64 + 49], arr_index[n4*64 + 53], n4%2==0);
//         compareSwap(arr[n4*64 + 50], arr[n4*64 + 54], arr_index[n4*64 + 50], arr_index[n4*64 + 54], n4%2==0);
//         compareSwap(arr[n4*64 + 51], arr[n4*64 + 55], arr_index[n4*64 + 51], arr_index[n4*64 + 55], n4%2==0);

//         compareSwap(arr[n4*64 + 56], arr[n4*64 + 60], arr_index[n4*64 + 56], arr_index[n4*64 + 60], n4%2==0);
//         compareSwap(arr[n4*64 + 57], arr[n4*64 + 61], arr_index[n4*64 + 57], arr_index[n4*64 + 61], n4%2==0);
//         compareSwap(arr[n4*64 + 58], arr[n4*64 + 62], arr_index[n4*64 + 58], arr_index[n4*64 + 62], n4%2==0);
//         compareSwap(arr[n4*64 + 59], arr[n4*64 + 63], arr_index[n4*64 + 59], arr_index[n4*64 + 63], n4%2==0);
//     }

//     for(int n5 = 0; n5 < 1; n5++){
// #pragma HLS UNROLL
//         compareSwap(arr[n5*64], arr[n5*64 + 2], arr_index[n5*64], arr_index[n5*64 + 2], n5%2==0);
//         compareSwap(arr[n5*64 + 1], arr[n5*64 + 3], arr_index[n5*64 + 1], arr_index[n5*64 + 3], n5%2==0);
//         compareSwap(arr[n5*64 + 4], arr[n5*64 + 6], arr_index[n5*64 + 4], arr_index[n5*64 + 6], n5%2==0);
//         compareSwap(arr[n5*64 + 5], arr[n5*64 + 7], arr_index[n5*64 + 5], arr_index[n5*64 + 7], n5%2==0);
//         compareSwap(arr[n5*64 + 8], arr[n5*64 + 10], arr_index[n5*64 + 8], arr_index[n5*64 + 10], n5%2==0);
//         compareSwap(arr[n5*64 + 9], arr[n5*64 + 11], arr_index[n5*64 + 9], arr_index[n5*64 + 11], n5%2==0);
//         compareSwap(arr[n5*64 + 12], arr[n5*64 + 14], arr_index[n5*64 + 12], arr_index[n5*64 + 14], n5%2==0);
//         compareSwap(arr[n5*64 + 13], arr[n5*64 + 15], arr_index[n5*64 + 13], arr_index[n5*64 + 15], n5%2==0);

//         compareSwap(arr[n5*64 + 16], arr[n5*64 + 18], arr_index[n5*64 + 16], arr_index[n5*64 + 18], n5%2==0);
//         compareSwap(arr[n5*64 + 17], arr[n5*64 + 19], arr_index[n5*64 + 17], arr_index[n5*64 + 19], n5%2==0);
//         compareSwap(arr[n5*64 + 20], arr[n5*64 + 22], arr_index[n5*64 + 20], arr_index[n5*64 + 22], n5%2==0);
//         compareSwap(arr[n5*64 + 21], arr[n5*64 + 23], arr_index[n5*64 + 21], arr_index[n5*64 + 23], n5%2==0);
//         compareSwap(arr[n5*64 + 24], arr[n5*64 + 26], arr_index[n5*64 + 24], arr_index[n5*64 + 26], n5%2==0);
//         compareSwap(arr[n5*64 + 25], arr[n5*64 + 27], arr_index[n5*64 + 25], arr_index[n5*64 + 27], n5%2==0);
//         compareSwap(arr[n5*64 + 28], arr[n5*64 + 30], arr_index[n5*64 + 28], arr_index[n5*64 + 30], n5%2==0);
//         compareSwap(arr[n5*64 + 29], arr[n5*64 + 31], arr_index[n5*64 + 29], arr_index[n5*64 + 31], n5%2==0);
//         compareSwap(arr[n5*64 + 32], arr[n5*64 + 34], arr_index[n5*64 + 32], arr_index[n5*64 + 34], n5%2==0);
//         compareSwap(arr[n5*64 + 33], arr[n5*64 + 35], arr_index[n5*64 + 33], arr_index[n5*64 + 35], n5%2==0);
        
//         compareSwap(arr[n5*64 + 36], arr[n5*64 + 38], arr_index[n5*64 + 36], arr_index[n5*64 + 38], n5%2==0);
//         compareSwap(arr[n5*64 + 37], arr[n5*64 + 39], arr_index[n5*64 + 37], arr_index[n5*64 + 39], n5%2==0);
//         compareSwap(arr[n5*64 + 40], arr[n5*64 + 42], arr_index[n5*64 + 40], arr_index[n5*64 + 42], n5%2==0);
//         compareSwap(arr[n5*64 + 41], arr[n5*64 + 43], arr_index[n5*64 + 41], arr_index[n5*64 + 43], n5%2==0);
//         compareSwap(arr[n5*64 + 44], arr[n5*64 + 46], arr_index[n5*64 + 44], arr_index[n5*64 + 46], n5%2==0);
//         compareSwap(arr[n5*64 + 45], arr[n5*64 + 47], arr_index[n5*64 + 45], arr_index[n5*64 + 47], n5%2==0);
//         compareSwap(arr[n5*64 + 48], arr[n5*64 + 50], arr_index[n5*64 + 48], arr_index[n5*64 + 50], n5%2==0);
//         compareSwap(arr[n5*64 + 49], arr[n5*64 + 51], arr_index[n5*64 + 49], arr_index[n5*64 + 51], n5%2==0);
//         compareSwap(arr[n5*64 + 52], arr[n5*64 + 54], arr_index[n5*64 + 52], arr_index[n5*64 + 54], n5%2==0);
//         compareSwap(arr[n5*64 + 53], arr[n5*64 + 55], arr_index[n5*64 + 53], arr_index[n5*64 + 55], n5%2==0);
//         compareSwap(arr[n5*64 + 56], arr[n5*64 + 58], arr_index[n5*64 + 56], arr_index[n5*64 + 58], n5%2==0);
//         compareSwap(arr[n5*64 + 57], arr[n5*64 + 59], arr_index[n5*64 + 57], arr_index[n5*64 + 59], n5%2==0);
//         compareSwap(arr[n5*64 + 60], arr[n5*64 + 62], arr_index[n5*64 + 60], arr_index[n5*64 + 62], n5%2==0);
//         compareSwap(arr[n5*64 + 61], arr[n5*64 + 63], arr_index[n5*64 + 61], arr_index[n5*64 + 63], n5%2==0);
//     }

//     for(int n6 = 0; n6 < 1; n6++){
// #pragma HLS UNROLL
//         compareSwap(arr[n6*64], arr[n6*64 + 1], arr_index[n6*64], arr_index[n6*64 + 1], n6%2==0);
//         compareSwap(arr[n6*64 + 2], arr[n6*64 + 3], arr_index[n6*64 + 2], arr_index[n6*64 + 3], n6%2==0);
//         compareSwap(arr[n6*64 + 4], arr[n6*64 + 5], arr_index[n6*64 + 4], arr_index[n6*64 + 5], n6%2==0);
//         compareSwap(arr[n6*64 + 6], arr[n6*64 + 7], arr_index[n6*64 + 6], arr_index[n6*64 + 7], n6%2==0);
//         compareSwap(arr[n6*64 + 8], arr[n6*64 + 9], arr_index[n6*64 + 8], arr_index[n6*64 + 9], n6%2==0);
//         compareSwap(arr[n6*64 + 10], arr[n6*64 + 11], arr_index[n6*64 + 10], arr_index[n6*64 + 11], n6%2==0);
//         compareSwap(arr[n6*64 + 12], arr[n6*64 + 13], arr_index[n6*64 + 12], arr_index[n6*64 + 13], n6%2==0);
//         compareSwap(arr[n6*64 + 14], arr[n6*64 + 15], arr_index[n6*64 + 14], arr_index[n6*64 + 15], n6%2==0);
//         compareSwap(arr[n6*64 + 16], arr[n6*64 + 17], arr_index[n6*64 + 16], arr_index[n6*64 + 17], n6%2==0);
//         compareSwap(arr[n6*64 + 18], arr[n6*64 + 19], arr_index[n6*64 + 18], arr_index[n6*64 + 19], n6%2==0);
//         compareSwap(arr[n6*64 + 20], arr[n6*64 + 21], arr_index[n6*64 + 20], arr_index[n6*64 + 21], n6%2==0);
//         compareSwap(arr[n6*64 + 22], arr[n6*64 + 23], arr_index[n6*64 + 22], arr_index[n6*64 + 23], n6%2==0);
//         compareSwap(arr[n6*64 + 24], arr[n6*64 + 25], arr_index[n6*64 + 24], arr_index[n6*64 + 25], n6%2==0);
//         compareSwap(arr[n6*64 + 26], arr[n6*64 + 27], arr_index[n6*64 + 26], arr_index[n6*64 + 27], n6%2==0);
//         compareSwap(arr[n6*64 + 28], arr[n6*64 + 29], arr_index[n6*64 + 28], arr_index[n6*64 + 29], n6%2==0);
//         compareSwap(arr[n6*64 + 30], arr[n6*64 + 31], arr_index[n6*64 + 30], arr_index[n6*64 + 31], n6%2==0);
//         compareSwap(arr[n6*64 + 32], arr[n6*64 + 33], arr_index[n6*64 + 32], arr_index[n6*64 + 33], n6%2==0);
//         compareSwap(arr[n6*64 + 34], arr[n6*64 + 35], arr_index[n6*64 + 34], arr_index[n6*64 + 35], n6%2==0);
//         compareSwap(arr[n6*64 + 36], arr[n6*64 + 37], arr_index[n6*64 + 36], arr_index[n6*64 + 37], n6%2==0);
//         compareSwap(arr[n6*64 + 38], arr[n6*64 + 39], arr_index[n6*64 + 38], arr_index[n6*64 + 39], n6%2==0);
//         compareSwap(arr[n6*64 + 40], arr[n6*64 + 41], arr_index[n6*64 + 40], arr_index[n6*64 + 41], n6%2==0);
//         compareSwap(arr[n6*64 + 42], arr[n6*64 + 43], arr_index[n6*64 + 42], arr_index[n6*64 + 43], n6%2==0);
//         compareSwap(arr[n6*64 + 44], arr[n6*64 + 45], arr_index[n6*64 + 44], arr_index[n6*64 + 45], n6%2==0);
//         compareSwap(arr[n6*64 + 46], arr[n6*64 + 47], arr_index[n6*64 + 46], arr_index[n6*64 + 47], n6%2==0);
//         compareSwap(arr[n6*64 + 48], arr[n6*64 + 49], arr_index[n6*64 + 48], arr_index[n6*64 + 49], n6%2==0);
//         compareSwap(arr[n6*64 + 50], arr[n6*64 + 51], arr_index[n6*64 + 50], arr_index[n6*64 + 51], n6%2==0);
//         compareSwap(arr[n6*64 + 52], arr[n6*64 + 53], arr_index[n6*64 + 52], arr_index[n6*64 + 53], n6%2==0);
//         compareSwap(arr[n6*64 + 54], arr[n6*64 + 55], arr_index[n6*64 + 54], arr_index[n6*64 + 55], n6%2==0);
//         compareSwap(arr[n6*64 + 56], arr[n6*64 + 57], arr_index[n6*64 + 56], arr_index[n6*64 + 57], n6%2==0);
//         compareSwap(arr[n6*64 + 58], arr[n6*64 + 59], arr_index[n6*64 + 58], arr_index[n6*64 + 59], n6%2==0);
//         compareSwap(arr[n6*64 + 60], arr[n6*64 + 61], arr_index[n6*64 + 60], arr_index[n6*64 + 61], n6%2==0);
//         compareSwap(arr[n6*64 + 62], arr[n6*64 + 63], arr_index[n6*64 + 62], arr_index[n6*64 + 63], n6%2==0);
//     }
    
}

void bitonicSort256_hw(type_dist_hw* arr, voxel_int* arr_index) 
{
// 对数组进行完全分区以实现并行访问
#pragma HLS ARRAY_PARTITION variable=arr complete
#pragma HLS ARRAY_PARTITION variable=arr_index complete
    //Stage 1
    for(int i = 0; i < 128; i++) {
#pragma HLS UNROLL factor = 16
        compareSwap(arr[i*2], arr[i*2 + 1], arr_index[i*2], arr_index[i*2 + 1], i%2==0);
    }

    //Stage 2
    for(int j = 0; j < 64; j++) {
#pragma HLS UNROLL factor = 8
        compareSwap(arr[j*4], arr[j*4 + 2], arr_index[j*4], arr_index[j*4 + 2], j%2==0);
        compareSwap(arr[j*4 + 1], arr[j*4 + 3], arr_index[j*4 + 1], arr_index[j*4 + 3], j%2==0);
    }
    for(int jj = 0; jj < 64; jj++) {
#pragma HLS UNROLL factor = 8
        compareSwap(arr[jj*4], arr[jj*4 + 1], arr_index[jj*4], arr_index[jj*4 + 1], jj%2==0);
        compareSwap(arr[jj*4 + 2], arr[jj*4 + 3], arr_index[jj*4 + 2], arr_index[jj*4 + 3], jj%2==0);
    }

    //Stage 3
    for(int k = 0; k < 32; k++) {
#pragma HLS UNROLL factor = 4
        compareSwap(arr[k*8], arr[k*8 + 4], arr_index[k*8], arr_index[k*8 + 4], k%2==0);
        compareSwap(arr[k*8 + 1], arr[k*8 + 5], arr_index[k*8 + 1], arr_index[k*8 + 5], k%2==0);
        compareSwap(arr[k*8 + 2], arr[k*8 + 6], arr_index[k*8 + 2], arr_index[k*8 + 6], k%2==0);
        compareSwap(arr[k*8 + 3], arr[k*8 + 7], arr_index[k*8 + 3], arr_index[k*8 + 7], k%2==0);
    }
    for(int kk = 0; kk < 32; kk++) {
#pragma HLS UNROLL factor = 4
        compareSwap(arr[kk*8], arr[kk*8 + 2], arr_index[kk*8], arr_index[kk*8 + 2], kk%2==0);
        compareSwap(arr[kk*8 + 1], arr[kk*8 + 3], arr_index[kk*8 + 1], arr_index[kk*8 + 3], kk%2==0);
        compareSwap(arr[kk*8 + 4], arr[kk*8 + 6], arr_index[kk*8 + 4], arr_index[kk*8 + 6], kk%2==0);
        compareSwap(arr[kk*8 + 5], arr[kk*8 + 7], arr_index[kk*8 + 5], arr_index[kk*8 + 7], kk%2==0);
    }
    for(int kkk = 0; kkk < 32; kkk++) {
#pragma HLS UNROLL factor = 4
        compareSwap(arr[kkk*8], arr[kkk*8 + 1], arr_index[kkk*8], arr_index[kkk*8 + 1], kkk%2==0);
        compareSwap(arr[kkk*8 + 2], arr[kkk*8 + 3], arr_index[kkk*8 + 2], arr_index[kkk*8 + 3], kkk%2==0);
        compareSwap(arr[kkk*8 + 4], arr[kkk*8 + 5], arr_index[kkk*8 + 4], arr_index[kkk*8 + 5], kkk%2==0);
        compareSwap(arr[kkk*8 + 6], arr[kkk*8 + 7], arr_index[kkk*8 + 6], arr_index[kkk*8 + 7], kkk%2==0);
    }

    //Stage 4
    for(int l = 0; l < 16; l++){
#pragma HLS UNROLL factor = 2
        compareSwap(arr[l*16], arr[l*16 + 8], arr_index[l*16], arr_index[l*16 + 8], l%2==0);
        compareSwap(arr[l*16 + 1], arr[l*16 + 9], arr_index[l*16 + 1], arr_index[l*16 + 9], l%2==0);
        compareSwap(arr[l*16 + 2], arr[l*16 + 10], arr_index[l*16 + 2], arr_index[l*16 + 10], l%2==0);
        compareSwap(arr[l*16 + 3], arr[l*16 + 11], arr_index[l*16 + 3], arr_index[l*16 + 11], l%2==0);
        compareSwap(arr[l*16 + 4], arr[l*16 + 12], arr_index[l*16 + 4], arr_index[l*16 + 12], l%2==0);
        compareSwap(arr[l*16 + 5], arr[l*16 + 13], arr_index[l*16 + 5], arr_index[l*16 + 13], l%2==0);
        compareSwap(arr[l*16 + 6], arr[l*16 + 14], arr_index[l*16 + 6], arr_index[l*16 + 14], l%2==0);
        compareSwap(arr[l*16 + 7], arr[l*16 + 15], arr_index[l*16 + 7], arr_index[l*16 + 15], l%2==0);
    }

    for(int ll = 0; ll < 16; ll++){
#pragma HLS UNROLL factor = 2
        compareSwap(arr[ll*16], arr[ll*16 + 4], arr_index[ll*16], arr_index[ll*16 + 4], ll%2==0);
        compareSwap(arr[ll*16 + 1], arr[ll*16 + 5], arr_index[ll*16 + 1], arr_index[ll*16 + 5], ll%2==0);
        compareSwap(arr[ll*16 + 2], arr[ll*16 + 6], arr_index[ll*16 + 2], arr_index[ll*16 + 6], ll%2==0);
        compareSwap(arr[ll*16 + 3], arr[ll*16 + 7], arr_index[ll*16 + 3], arr_index[ll*16 + 7], ll%2==0);
        compareSwap(arr[ll*16 + 8], arr[ll*16 + 12], arr_index[ll*16 + 8], arr_index[ll*16 + 12], ll%2==0);
        compareSwap(arr[ll*16 + 9], arr[ll*16 + 13], arr_index[ll*16 + 9], arr_index[ll*16 + 13], ll%2==0);
        compareSwap(arr[ll*16 + 10], arr[ll*16 + 14], arr_index[ll*16 + 10], arr_index[ll*16 + 14], ll%2==0);
        compareSwap(arr[ll*16 + 11], arr[ll*16 + 15], arr_index[ll*16 + 11], arr_index[ll*16 + 15], ll%2==0);
    }

    for (int lll = 0; lll < 16; lll++){
#pragma HLS UNROLL factor = 2
        compareSwap(arr[lll*16], arr[lll*16 + 2], arr_index[lll*16], arr_index[lll*16 + 2], lll%2==0);
        compareSwap(arr[lll*16 + 1], arr[lll*16 + 3], arr_index[lll*16 + 1], arr_index[lll*16 + 3], lll%2==0);
        compareSwap(arr[lll*16 + 4], arr[lll*16 + 6], arr_index[lll*16 + 4], arr_index[lll*16 + 6], lll%2==0);
        compareSwap(arr[lll*16 + 5], arr[lll*16 + 7], arr_index[lll*16 + 5], arr_index[lll*16 + 7], lll%2==0);
        compareSwap(arr[lll*16 + 8], arr[lll*16 + 10], arr_index[lll*16 + 8], arr_index[lll*16 + 10], lll%2==0);
        compareSwap(arr[lll*16 + 9], arr[lll*16 + 11], arr_index[lll*16 + 9], arr_index[lll*16 + 11], lll%2==0);
        compareSwap(arr[lll*16 + 12], arr[lll*16 + 14], arr_index[lll*16 + 12], arr_index[lll*16 + 14], lll%2==0);
        compareSwap(arr[lll*16 + 13], arr[lll*16 + 15], arr_index[lll*16 + 13], arr_index[lll*16 + 15], lll%2==0);
    }

    for (int llll = 0; llll < 16; llll++){
#pragma HLS UNROLL factor = 2
        compareSwap(arr[llll*16], arr[llll*16 + 1], arr_index[llll*16], arr_index[llll*16 + 1], llll%2==0);
        compareSwap(arr[llll*16 + 2], arr[llll*16 + 3], arr_index[llll*16 + 2], arr_index[llll*16 + 3], llll%2==0);
        compareSwap(arr[llll*16 + 4], arr[llll*16 + 5], arr_index[llll*16 + 4], arr_index[llll*16 + 5], llll%2==0);
        compareSwap(arr[llll*16 + 6], arr[llll*16 + 7], arr_index[llll*16 + 6], arr_index[llll*16 + 7], llll%2==0);
        compareSwap(arr[llll*16 + 8], arr[llll*16 + 9], arr_index[llll*16 + 8], arr_index[llll*16 + 9], llll%2==0);
        compareSwap(arr[llll*16 + 10], arr[llll*16 + 11], arr_index[llll*16 + 10], arr_index[llll*16 + 11], llll%2==0);
        compareSwap(arr[llll*16 + 12], arr[llll*16 + 13], arr_index[llll*16 + 12], arr_index[llll*16 + 13], llll%2==0);
        compareSwap(arr[llll*16 + 14], arr[llll*16 + 15], arr_index[llll*16 + 14], arr_index[llll*16 + 15], llll%2==0);
    }

    //Stage 5

    for(int m = 0; m < 8; m++){
// #pragma HLS UNROLL
        compareSwap(arr[m*32], arr[m*32 + 16], arr_index[m*32], arr_index[m*32 + 16], m%2==0);
        compareSwap(arr[m*32 + 1], arr[m*32 + 17], arr_index[m*32 + 1], arr_index[m*32 + 17], m%2==0);
        compareSwap(arr[m*32 + 2], arr[m*32 + 18], arr_index[m*32 + 2], arr_index[m*32 + 18], m%2==0);
        compareSwap(arr[m*32 + 3], arr[m*32 + 19], arr_index[m*32 + 3], arr_index[m*32 + 19], m%2==0);
        compareSwap(arr[m*32 + 4], arr[m*32 + 20], arr_index[m*32 + 4], arr_index[m*32 + 20], m%2==0);
        compareSwap(arr[m*32 + 5], arr[m*32 + 21], arr_index[m*32 + 5], arr_index[m*32 + 21], m%2==0);
        compareSwap(arr[m*32 + 6], arr[m*32 + 22], arr_index[m*32 + 6], arr_index[m*32 + 22], m%2==0);
        compareSwap(arr[m*32 + 7], arr[m*32 + 23], arr_index[m*32 + 7], arr_index[m*32 + 23], m%2==0);
        compareSwap(arr[m*32 + 8], arr[m*32 + 24], arr_index[m*32 + 8], arr_index[m*32 + 24], m%2==0);
        compareSwap(arr[m*32 + 9], arr[m*32 + 25], arr_index[m*32 + 9], arr_index[m*32 + 25], m%2==0);
        compareSwap(arr[m*32 + 10], arr[m*32 + 26], arr_index[m*32 + 10], arr_index[m*32 + 26], m%2==0);
        compareSwap(arr[m*32 + 11], arr[m*32 + 27], arr_index[m*32 + 11], arr_index[m*32 + 27], m%2==0);
        compareSwap(arr[m*32 + 12], arr[m*32 + 28], arr_index[m*32 + 12], arr_index[m*32 + 28], m%2==0);
        compareSwap(arr[m*32 + 13], arr[m*32 + 29], arr_index[m*32 + 13], arr_index[m*32 + 29], m%2==0);
        compareSwap(arr[m*32 + 14], arr[m*32 + 30], arr_index[m*32 + 14], arr_index[m*32 + 30], m%2==0);
        compareSwap(arr[m*32 + 15], arr[m*32 + 31], arr_index[m*32 + 15], arr_index[m*32 + 31], m%2==0);
    }

    for(int mm = 0; mm < 8; mm++){
// #pragma HLS UNROLL
        compareSwap(arr[mm*32], arr[mm*32 + 8], arr_index[mm*32], arr_index[mm*32 + 8], mm%2==0);
        compareSwap(arr[mm*32 + 1], arr[mm*32 + 9], arr_index[mm*32 + 1], arr_index[mm*32 + 9], mm%2==0);
        compareSwap(arr[mm*32 + 2], arr[mm*32 + 10], arr_index[mm*32 + 2], arr_index[mm*32 + 10], mm%2==0);
        compareSwap(arr[mm*32 + 3], arr[mm*32 + 11], arr_index[mm*32 + 3], arr_index[mm*32 + 11], mm%2==0);
        compareSwap(arr[mm*32 + 4], arr[mm*32 + 12], arr_index[mm*32 + 4], arr_index[mm*32 + 12], mm%2==0);
        compareSwap(arr[mm*32 + 5], arr[mm*32 + 13], arr_index[mm*32 + 5], arr_index[mm*32 + 13], mm%2==0);
        compareSwap(arr[mm*32 + 6], arr[mm*32 + 14], arr_index[mm*32 + 6], arr_index[mm*32 + 14], mm%2==0);
        compareSwap(arr[mm*32 + 7], arr[mm*32 + 15], arr_index[mm*32 + 7], arr_index[mm*32 + 15], mm%2==0);
        compareSwap(arr[mm*32 + 16], arr[mm*32 + 24], arr_index[mm*32 + 16], arr_index[mm*32 + 24], mm%2==0);
        compareSwap(arr[mm*32 + 17], arr[mm*32 + 25], arr_index[mm*32 + 17], arr_index[mm*32 + 25], mm%2==0);
        compareSwap(arr[mm*32 + 18], arr[mm*32 + 26], arr_index[mm*32 + 18], arr_index[mm*32 + 26], mm%2==0);
        compareSwap(arr[mm*32 + 19], arr[mm*32 + 27], arr_index[mm*32 + 19], arr_index[mm*32 + 27], mm%2==0);
        compareSwap(arr[mm*32 + 20], arr[mm*32 + 28], arr_index[mm*32 + 20], arr_index[mm*32 + 28], mm%2==0);
        compareSwap(arr[mm*32 + 21], arr[mm*32 + 29], arr_index[mm*32 + 21], arr_index[mm*32 + 29], mm%2==0);
        compareSwap(arr[mm*32 + 22], arr[mm*32 + 30], arr_index[mm*32 + 22], arr_index[mm*32 + 30], mm%2==0);
        compareSwap(arr[mm*32 + 23], arr[mm*32 + 31], arr_index[mm*32 + 23], arr_index[mm*32 + 31], mm%2==0);
    }

    for(int mmm = 0; mmm < 8; mmm++){
// #pragma HLS UNROLL
        compareSwap(arr[mmm*32], arr[mmm*32 + 4], arr_index[mmm*32], arr_index[mmm*32 + 4], mmm%2==0);
        compareSwap(arr[mmm*32 + 1], arr[mmm*32 + 5], arr_index[mmm*32 + 1], arr_index[mmm*32 + 5], mmm%2==0);
        compareSwap(arr[mmm*32 + 2], arr[mmm*32 + 6], arr_index[mmm*32 + 2], arr_index[mmm*32 + 6], mmm%2==0);
        compareSwap(arr[mmm*32 + 3], arr[mmm*32 + 7], arr_index[mmm*32 + 3], arr_index[mmm*32 + 7], mmm%2==0);
        compareSwap(arr[mmm*32 + 8], arr[mmm*32 + 12], arr_index[mmm*32 + 8], arr_index[mmm*32 + 12], mmm%2==0);
        compareSwap(arr[mmm*32 + 9], arr[mmm*32 + 13], arr_index[mmm*32 + 9], arr_index[mmm*32 + 13], mmm%2==0);
        compareSwap(arr[mmm*32 + 10], arr[mmm*32 + 14], arr_index[mmm*32 + 10], arr_index[mmm*32 + 14], mmm%2==0);
        compareSwap(arr[mmm*32 + 11], arr[mmm*32 + 15], arr_index[mmm*32 + 11], arr_index[mmm*32 + 15], mmm%2==0);
        compareSwap(arr[mmm*32 + 16], arr[mmm*32 + 20], arr_index[mmm*32 + 16], arr_index[mmm*32 + 20], mmm%2==0);
        compareSwap(arr[mmm*32 + 17], arr[mmm*32 + 21], arr_index[mmm*32 + 17], arr_index[mmm*32 + 21], mmm%2==0);
        compareSwap(arr[mmm*32 + 18], arr[mmm*32 + 22], arr_index[mmm*32 + 18], arr_index[mmm*32 + 22], mmm%2==0);
        compareSwap(arr[mmm*32 + 19], arr[mmm*32 + 23], arr_index[mmm*32 + 19], arr_index[mmm*32 + 23], mmm%2==0);
        compareSwap(arr[mmm*32 + 24], arr[mmm*32 + 28], arr_index[mmm*32 + 24], arr_index[mmm*32 + 28], mmm%2==0);
        compareSwap(arr[mmm*32 + 25], arr[mmm*32 + 29], arr_index[mmm*32 + 25], arr_index[mmm*32 + 29], mmm%2==0);
        compareSwap(arr[mmm*32 + 26], arr[mmm*32 + 30], arr_index[mmm*32 + 26], arr_index[mmm*32 + 30], mmm%2==0);
        compareSwap(arr[mmm*32 + 27], arr[mmm*32 + 31], arr_index[mmm*32 + 27], arr_index[mmm*32 + 31], mmm%2==0);
    }

    for(int mmmm = 0; mmmm < 8; mmmm ++){
// #pragma HLS UNROLL
        compareSwap(arr[mmmm*32], arr[mmmm*32 + 2], arr_index[mmmm*32], arr_index[mmmm*32 + 2], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 1], arr[mmmm*32 + 3], arr_index[mmmm*32 + 1], arr_index[mmmm*32 + 3], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 4], arr[mmmm*32 + 6], arr_index[mmmm*32 + 4], arr_index[mmmm*32 + 6], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 5], arr[mmmm*32 + 7], arr_index[mmmm*32 + 5], arr_index[mmmm*32 + 7], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 8], arr[mmmm*32 + 10], arr_index[mmmm*32 + 8], arr_index[mmmm*32 + 10], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 9], arr[mmmm*32 + 11], arr_index[mmmm*32 + 9], arr_index[mmmm*32 + 11], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 12], arr[mmmm*32 + 14], arr_index[mmmm*32 + 12], arr_index[mmmm*32 + 14], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 13], arr[mmmm*32 + 15], arr_index[mmmm*32 + 13], arr_index[mmmm*32 + 15], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 16], arr[mmmm*32 + 18], arr_index[mmmm*32 + 16], arr_index[mmmm*32 + 18], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 17], arr[mmmm*32 + 19], arr_index[mmmm*32 + 17], arr_index[mmmm*32 + 19], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 20], arr[mmmm*32 + 22], arr_index[mmmm*32 + 20], arr_index[mmmm*32 + 22], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 21], arr[mmmm*32 + 23], arr_index[mmmm*32 + 21], arr_index[mmmm*32 + 23], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 24], arr[mmmm*32 + 26], arr_index[mmmm*32 + 24], arr_index[mmmm*32 + 26], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 25], arr[mmmm*32 + 27], arr_index[mmmm*32 + 25], arr_index[mmmm*32 + 27], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 28], arr[mmmm*32 + 30], arr_index[mmmm*32 + 28], arr_index[mmmm*32 + 30], mmmm%2==0);
        compareSwap(arr[mmmm*32 + 29], arr[mmmm*32 + 31], arr_index[mmmm*32 + 29], arr_index[mmmm*32 + 31], mmmm%2==0);
    }

    for(int m5 = 0; m5 < 8; m5++){
// #pragma HLS UNROLL
        compareSwap(arr[m5*32], arr[m5*32 + 1], arr_index[m5*32], arr_index[m5*32 + 1], m5%2==0);
        compareSwap(arr[m5*32 + 2], arr[m5*32 + 3], arr_index[m5*32 + 2], arr_index[m5*32 + 3], m5%2==0);
        compareSwap(arr[m5*32 + 4], arr[m5*32 + 5], arr_index[m5*32 + 4], arr_index[m5*32 + 5], m5%2==0);
        compareSwap(arr[m5*32 + 6], arr[m5*32 + 7], arr_index[m5*32 + 6], arr_index[m5*32 + 7], m5%2==0);
        compareSwap(arr[m5*32 + 8], arr[m5*32 + 9], arr_index[m5*32 + 8], arr_index[m5*32 + 9], m5%2==0);
        compareSwap(arr[m5*32 + 10], arr[m5*32 + 11], arr_index[m5*32 + 10], arr_index[m5*32 + 11], m5%2==0);
        compareSwap(arr[m5*32 + 12], arr[m5*32 + 13], arr_index[m5*32 + 12], arr_index[m5*32 + 13], m5%2==0);
        compareSwap(arr[m5*32 + 14], arr[m5*32 + 15], arr_index[m5*32 + 14], arr_index[m5*32 + 15], m5%2==0);
        compareSwap(arr[m5*32 + 16], arr[m5*32 + 17], arr_index[m5*32 + 16], arr_index[m5*32 + 17], m5%2==0);
        compareSwap(arr[m5*32 + 18], arr[m5*32 + 19], arr_index[m5*32 + 18], arr_index[m5*32 + 19], m5%2==0);
        compareSwap(arr[m5*32 + 20], arr[m5*32 + 21], arr_index[m5*32 + 20], arr_index[m5*32 + 21], m5%2==0);
        compareSwap(arr[m5*32 + 22], arr[m5*32 + 23], arr_index[m5*32 + 22], arr_index[m5*32 + 23], m5%2==0);
        compareSwap(arr[m5*32 + 24], arr[m5*32 + 25], arr_index[m5*32 + 24], arr_index[m5*32 + 25], m5%2==0);
        compareSwap(arr[m5*32 + 26], arr[m5*32 + 27], arr_index[m5*32 + 26], arr_index[m5*32 + 27], m5%2==0);
        compareSwap(arr[m5*32 + 28], arr[m5*32 + 29], arr_index[m5*32 + 28], arr_index[m5*32 + 29], m5%2==0);
        compareSwap(arr[m5*32 + 30], arr[m5*32 + 31], arr_index[m5*32 + 30], arr_index[m5*32 + 31], m5%2==0);
    }

    //Stage 6

    for(int n = 0; n < 4; n++){

        compareSwap(arr[n*64], arr[n*64 + 32], arr_index[n*64], arr_index[n*64 + 32], n%2==0);
        compareSwap(arr[n*64 + 1], arr[n*64 + 33], arr_index[n*64 + 1], arr_index[n*64 + 33], n%2==0);
        compareSwap(arr[n*64 + 2], arr[n*64 + 34], arr_index[n*64 + 2], arr_index[n*64 + 34], n%2==0);
        compareSwap(arr[n*64 + 3], arr[n*64 + 35], arr_index[n*64 + 3], arr_index[n*64 + 35], n%2==0);
        compareSwap(arr[n*64 + 4], arr[n*64 + 36], arr_index[n*64 + 4], arr_index[n*64 + 36], n%2==0);
        compareSwap(arr[n*64 + 5], arr[n*64 + 37], arr_index[n*64 + 5], arr_index[n*64 + 37], n%2==0);
        compareSwap(arr[n*64 + 6], arr[n*64 + 38], arr_index[n*64 + 6], arr_index[n*64 + 38], n%2==0);
        compareSwap(arr[n*64 + 7], arr[n*64 + 39], arr_index[n*64 + 7], arr_index[n*64 + 39], n%2==0);
        compareSwap(arr[n*64 + 8], arr[n*64 + 40], arr_index[n*64 + 8], arr_index[n*64 + 40], n%2==0);
        compareSwap(arr[n*64 + 9], arr[n*64 + 41], arr_index[n*64 + 9], arr_index[n*64 + 41], n%2==0);
        compareSwap(arr[n*64 + 10], arr[n*64 + 42], arr_index[n*64 + 10], arr_index[n*64 + 42], n%2==0);
        compareSwap(arr[n*64 + 11], arr[n*64 + 43], arr_index[n*64 + 11], arr_index[n*64 + 43], n%2==0);
        compareSwap(arr[n*64 + 12], arr[n*64 + 44], arr_index[n*64 + 12], arr_index[n*64 + 44], n%2==0);
        compareSwap(arr[n*64 + 13], arr[n*64 + 45], arr_index[n*64 + 13], arr_index[n*64 + 45], n%2==0);
        compareSwap(arr[n*64 + 14], arr[n*64 + 46], arr_index[n*64 + 14], arr_index[n*64 + 46], n%2==0);
        compareSwap(arr[n*64 + 15], arr[n*64 + 47], arr_index[n*64 + 15], arr_index[n*64 + 47], n%2==0);
        compareSwap(arr[n*64 + 16], arr[n*64 + 48], arr_index[n*64 + 16], arr_index[n*64 + 48], n%2==0);
        compareSwap(arr[n*64 + 17], arr[n*64 + 49], arr_index[n*64 + 17], arr_index[n*64 + 49], n%2==0);
        compareSwap(arr[n*64 + 18], arr[n*64 + 50], arr_index[n*64 + 18], arr_index[n*64 + 50], n%2==0);
        compareSwap(arr[n*64 + 19], arr[n*64 + 51], arr_index[n*64 + 19], arr_index[n*64 + 51], n%2==0);
        compareSwap(arr[n*64 + 20], arr[n*64 + 52], arr_index[n*64 + 20], arr_index[n*64 + 52], n%2==0);
        compareSwap(arr[n*64 + 21], arr[n*64 + 53], arr_index[n*64 + 21], arr_index[n*64 + 53], n%2==0);
        compareSwap(arr[n*64 + 22], arr[n*64 + 54], arr_index[n*64 + 22], arr_index[n*64 + 54], n%2==0);
        compareSwap(arr[n*64 + 23], arr[n*64 + 55], arr_index[n*64 + 23], arr_index[n*64 + 55], n%2==0);
        compareSwap(arr[n*64 + 24], arr[n*64 + 56], arr_index[n*64 + 24], arr_index[n*64 + 56], n%2==0);
        compareSwap(arr[n*64 + 25], arr[n*64 + 57], arr_index[n*64 + 25], arr_index[n*64 + 57], n%2==0);
        compareSwap(arr[n*64 + 26], arr[n*64 + 58], arr_index[n*64 + 26], arr_index[n*64 + 58], n%2==0);
        compareSwap(arr[n*64 + 27], arr[n*64 + 59], arr_index[n*64 + 27], arr_index[n*64 + 59], n%2==0);
        compareSwap(arr[n*64 + 28], arr[n*64 + 60], arr_index[n*64 + 28], arr_index[n*64 + 60], n%2==0);
        compareSwap(arr[n*64 + 29], arr[n*64 + 61], arr_index[n*64 + 29], arr_index[n*64 + 61], n%2==0);
        compareSwap(arr[n*64 + 30], arr[n*64 + 62], arr_index[n*64 + 30], arr_index[n*64 + 62], n%2==0);
        compareSwap(arr[n*64 + 31], arr[n*64 + 63], arr_index[n*64 + 31], arr_index[n*64 + 63], n%2==0);
    }

    for(int n2 = 0; n2 < 4; n2++){

        compareSwap(arr[n2*64], arr[n2*64 + 16], arr_index[n2*64], arr_index[n2*64 + 16], n2%2==0);
        compareSwap(arr[n2*64 + 1], arr[n2*64 + 17], arr_index[n2*64 + 1], arr_index[n2*64 + 17], n2%2==0);
        compareSwap(arr[n2*64 + 2], arr[n2*64 + 18], arr_index[n2*64 + 2], arr_index[n2*64 + 18], n2%2==0);
        compareSwap(arr[n2*64 + 3], arr[n2*64 + 19], arr_index[n2*64 + 3], arr_index[n2*64 + 19], n2%2==0);
        compareSwap(arr[n2*64 + 4], arr[n2*64 + 20], arr_index[n2*64 + 4], arr_index[n2*64 + 20], n2%2==0);
        compareSwap(arr[n2*64 + 5], arr[n2*64 + 21], arr_index[n2*64 + 5], arr_index[n2*64 + 21], n2%2==0);
        compareSwap(arr[n2*64 + 6], arr[n2*64 + 22], arr_index[n2*64 + 6], arr_index[n2*64 + 22], n2%2==0);
        compareSwap(arr[n2*64 + 7], arr[n2*64 + 23], arr_index[n2*64 + 7], arr_index[n2*64 + 23], n2%2==0);
        compareSwap(arr[n2*64 + 8], arr[n2*64 + 24], arr_index[n2*64 + 8], arr_index[n2*64 + 24], n2%2==0);
        compareSwap(arr[n2*64 + 9], arr[n2*64 + 25], arr_index[n2*64 + 9], arr_index[n2*64 + 25], n2%2==0);
        compareSwap(arr[n2*64 + 10], arr[n2*64 + 26], arr_index[n2*64 + 10], arr_index[n2*64 + 26], n2%2==0);
        compareSwap(arr[n2*64 + 11], arr[n2*64 + 27], arr_index[n2*64 + 11], arr_index[n2*64 + 27], n2%2==0);
        compareSwap(arr[n2*64 + 12], arr[n2*64 + 28], arr_index[n2*64 + 12], arr_index[n2*64 + 28], n2%2==0);
        compareSwap(arr[n2*64 + 13], arr[n2*64 + 29], arr_index[n2*64 + 13], arr_index[n2*64 + 29], n2%2==0);
        compareSwap(arr[n2*64 + 14], arr[n2*64 + 30], arr_index[n2*64 + 14], arr_index[n2*64 + 30], n2%2==0);
        compareSwap(arr[n2*64 + 15], arr[n2*64 + 31], arr_index[n2*64 + 15], arr_index[n2*64 + 31], n2%2==0);
        
        compareSwap(arr[n2*64 + 32], arr[n2*64 + 48], arr_index[n2*64 + 32], arr_index[n2*64 + 48], n2%2==0);
        compareSwap(arr[n2*64 + 33], arr[n2*64 + 49], arr_index[n2*64 + 33], arr_index[n2*64 + 49], n2%2==0);
        compareSwap(arr[n2*64 + 34], arr[n2*64 + 50], arr_index[n2*64 + 34], arr_index[n2*64 + 50], n2%2==0);
        compareSwap(arr[n2*64 + 35], arr[n2*64 + 51], arr_index[n2*64 + 35], arr_index[n2*64 + 51], n2%2==0);
        compareSwap(arr[n2*64 + 36], arr[n2*64 + 52], arr_index[n2*64 + 36], arr_index[n2*64 + 52], n2%2==0);
        compareSwap(arr[n2*64 + 37], arr[n2*64 + 53], arr_index[n2*64 + 37], arr_index[n2*64 + 53], n2%2==0);
        compareSwap(arr[n2*64 + 38], arr[n2*64 + 54], arr_index[n2*64 + 38], arr_index[n2*64 + 54], n2%2==0);
        compareSwap(arr[n2*64 + 39], arr[n2*64 + 55], arr_index[n2*64 + 39], arr_index[n2*64 + 55], n2%2==0);
        compareSwap(arr[n2*64 + 40], arr[n2*64 + 56], arr_index[n2*64 + 40], arr_index[n2*64 + 56], n2%2==0);
        compareSwap(arr[n2*64 + 41], arr[n2*64 + 57], arr_index[n2*64 + 41], arr_index[n2*64 + 57], n2%2==0);
        compareSwap(arr[n2*64 + 42], arr[n2*64 + 58], arr_index[n2*64 + 42], arr_index[n2*64 + 58], n2%2==0);
        compareSwap(arr[n2*64 + 43], arr[n2*64 + 59], arr_index[n2*64 + 43], arr_index[n2*64 + 59], n2%2==0);
        compareSwap(arr[n2*64 + 44], arr[n2*64 + 60], arr_index[n2*64 + 44], arr_index[n2*64 + 60], n2%2==0);
        compareSwap(arr[n2*64 + 45], arr[n2*64 + 61], arr_index[n2*64 + 45], arr_index[n2*64 + 61], n2%2==0);
        compareSwap(arr[n2*64 + 46], arr[n2*64 + 62], arr_index[n2*64 + 46], arr_index[n2*64 + 62], n2%2==0);
        compareSwap(arr[n2*64 + 47], arr[n2*64 + 63], arr_index[n2*64 + 47], arr_index[n2*64 + 63], n2%2==0);
    }

    for(int n3 = 0; n3 < 4; n3++){

        compareSwap(arr[n3*64], arr[n3*64 + 8], arr_index[n3*64], arr_index[n3*64 + 8], n3%2==0);
        compareSwap(arr[n3*64 + 1], arr[n3*64 + 9], arr_index[n3*64 + 1], arr_index[n3*64 + 9], n3%2==0);
        compareSwap(arr[n3*64 + 2], arr[n3*64 + 10], arr_index[n3*64 + 2], arr_index[n3*64 + 10], n3%2==0);
        compareSwap(arr[n3*64 + 3], arr[n3*64 + 11], arr_index[n3*64 + 3], arr_index[n3*64 + 11], n3%2==0);
        compareSwap(arr[n3*64 + 4], arr[n3*64 + 12], arr_index[n3*64 + 4], arr_index[n3*64 + 12], n3%2==0);
        compareSwap(arr[n3*64 + 5], arr[n3*64 + 13], arr_index[n3*64 + 5], arr_index[n3*64 + 13], n3%2==0);
        compareSwap(arr[n3*64 + 6], arr[n3*64 + 14], arr_index[n3*64 + 6], arr_index[n3*64 + 14], n3%2==0);
        compareSwap(arr[n3*64 + 7], arr[n3*64 + 15], arr_index[n3*64 + 7], arr_index[n3*64 + 15], n3%2==0);

        compareSwap(arr[n3*64 + 16], arr[n3*64 + 24], arr_index[n3*64 + 16], arr_index[n3*64 + 24], n3%2==0);
        compareSwap(arr[n3*64 + 17], arr[n3*64 + 25], arr_index[n3*64 + 17], arr_index[n3*64 + 25], n3%2==0);
        compareSwap(arr[n3*64 + 18], arr[n3*64 + 26], arr_index[n3*64 + 18], arr_index[n3*64 + 26], n3%2==0);
        compareSwap(arr[n3*64 + 19], arr[n3*64 + 27], arr_index[n3*64 + 19], arr_index[n3*64 + 27], n3%2==0);
        compareSwap(arr[n3*64 + 20], arr[n3*64 + 28], arr_index[n3*64 + 20], arr_index[n3*64 + 28], n3%2==0);
        compareSwap(arr[n3*64 + 21], arr[n3*64 + 29], arr_index[n3*64 + 21], arr_index[n3*64 + 29], n3%2==0);
        compareSwap(arr[n3*64 + 22], arr[n3*64 + 30], arr_index[n3*64 + 22], arr_index[n3*64 + 30], n3%2==0);
        compareSwap(arr[n3*64 + 23], arr[n3*64 + 31], arr_index[n3*64 + 23], arr_index[n3*64 + 31], n3%2==0);

        compareSwap(arr[n3*64 + 32], arr[n3*64 + 40], arr_index[n3*64 + 32], arr_index[n3*64 + 40], n3%2==0);
        compareSwap(arr[n3*64 + 33], arr[n3*64 + 41], arr_index[n3*64 + 33], arr_index[n3*64 + 41], n3%2==0);
        compareSwap(arr[n3*64 + 34], arr[n3*64 + 42], arr_index[n3*64 + 34], arr_index[n3*64 + 42], n3%2==0);
        compareSwap(arr[n3*64 + 35], arr[n3*64 + 43], arr_index[n3*64 + 35], arr_index[n3*64 + 43], n3%2==0);
        compareSwap(arr[n3*64 + 36], arr[n3*64 + 44], arr_index[n3*64 + 36], arr_index[n3*64 + 44], n3%2==0);
        compareSwap(arr[n3*64 + 37], arr[n3*64 + 45], arr_index[n3*64 + 37], arr_index[n3*64 + 45], n3%2==0);
        compareSwap(arr[n3*64 + 38], arr[n3*64 + 46], arr_index[n3*64 + 38], arr_index[n3*64 + 46], n3%2==0);
        compareSwap(arr[n3*64 + 39], arr[n3*64 + 47], arr_index[n3*64 + 39], arr_index[n3*64 + 47], n3%2==0);
        
        compareSwap(arr[n3*64 + 48], arr[n3*64 + 56], arr_index[n3*64 + 48], arr_index[n3*64 + 56], n3%2==0);
        compareSwap(arr[n3*64 + 49], arr[n3*64 + 57], arr_index[n3*64 + 49], arr_index[n3*64 + 57], n3%2==0);
        compareSwap(arr[n3*64 + 50], arr[n3*64 + 58], arr_index[n3*64 + 50], arr_index[n3*64 + 58], n3%2==0);
        compareSwap(arr[n3*64 + 51], arr[n3*64 + 59], arr_index[n3*64 + 51], arr_index[n3*64 + 59], n3%2==0);
        compareSwap(arr[n3*64 + 52], arr[n3*64 + 60], arr_index[n3*64 + 52], arr_index[n3*64 + 60], n3%2==0);
        compareSwap(arr[n3*64 + 53], arr[n3*64 + 61], arr_index[n3*64 + 53], arr_index[n3*64 + 61], n3%2==0);
        compareSwap(arr[n3*64 + 54], arr[n3*64 + 62], arr_index[n3*64 + 54], arr_index[n3*64 + 62], n3%2==0);
        compareSwap(arr[n3*64 + 55], arr[n3*64 + 63], arr_index[n3*64 + 55], arr_index[n3*64 + 63], n3%2==0);
    }

    for(int n4 = 0; n4 < 4; n4++){

        compareSwap(arr[n4*64], arr[n4*64 + 4], arr_index[n4*64], arr_index[n4*64 + 4], n4%2==0);
        compareSwap(arr[n4*64 + 1], arr[n4*64 + 5], arr_index[n4*64 + 1], arr_index[n4*64 + 5], n4%2==0);
        compareSwap(arr[n4*64 + 2], arr[n4*64 + 6], arr_index[n4*64 + 2], arr_index[n4*64 + 6], n4%2==0);
        compareSwap(arr[n4*64 + 3], arr[n4*64 + 7], arr_index[n4*64 + 3], arr_index[n4*64 + 7], n4%2==0);
        compareSwap(arr[n4*64 + 8], arr[n4*64 + 12], arr_index[n4*64 + 8], arr_index[n4*64 + 12], n4%2==0);
        compareSwap(arr[n4*64 + 9], arr[n4*64 + 13], arr_index[n4*64 + 9], arr_index[n4*64 + 13], n4%2==0);
        compareSwap(arr[n4*64 + 10], arr[n4*64 + 14], arr_index[n4*64 + 10], arr_index[n4*64 + 14], n4%2==0);
        compareSwap(arr[n4*64 + 11], arr[n4*64 + 15], arr_index[n4*64 + 11], arr_index[n4*64 + 15], n4%2==0);
        compareSwap(arr[n4*64 + 16], arr[n4*64 + 20], arr_index[n4*64 + 16], arr_index[n4*64 + 20], n4%2==0);
        compareSwap(arr[n4*64 + 17], arr[n4*64 + 21], arr_index[n4*64 + 17], arr_index[n4*64 + 21], n4%2==0);
        compareSwap(arr[n4*64 + 18], arr[n4*64 + 22], arr_index[n4*64 + 18], arr_index[n4*64 + 22], n4%2==0);
        compareSwap(arr[n4*64 + 19], arr[n4*64 + 23], arr_index[n4*64 + 19], arr_index[n4*64 + 23], n4%2==0);
        compareSwap(arr[n4*64 + 24], arr[n4*64 + 28], arr_index[n4*64 + 24], arr_index[n4*64 + 28], n4%2==0);
        compareSwap(arr[n4*64 + 25], arr[n4*64 + 29], arr_index[n4*64 + 25], arr_index[n4*64 + 29], n4%2==0);
        compareSwap(arr[n4*64 + 26], arr[n4*64 + 30], arr_index[n4*64 + 26], arr_index[n4*64 + 30], n4%2==0);
        compareSwap(arr[n4*64 + 27], arr[n4*64 + 31], arr_index[n4*64 + 27], arr_index[n4*64 + 31], n4%2==0);
        compareSwap(arr[n4*64 + 32], arr[n4*64 + 36], arr_index[n4*64 + 32], arr_index[n4*64 + 36], n4%2==0);
        compareSwap(arr[n4*64 + 33], arr[n4*64 + 37], arr_index[n4*64 + 33], arr_index[n4*64 + 37], n4%2==0);
        compareSwap(arr[n4*64 + 34], arr[n4*64 + 38], arr_index[n4*64 + 34], arr_index[n4*64 + 38], n4%2==0);
        compareSwap(arr[n4*64 + 35], arr[n4*64 + 39], arr_index[n4*64 + 35], arr_index[n4*64 + 39], n4%2==0);

        compareSwap(arr[n4*64 + 40], arr[n4*64 + 44], arr_index[n4*64 + 40], arr_index[n4*64 + 44], n4%2==0);
        compareSwap(arr[n4*64 + 41], arr[n4*64 + 45], arr_index[n4*64 + 41], arr_index[n4*64 + 45], n4%2==0);
        compareSwap(arr[n4*64 + 42], arr[n4*64 + 46], arr_index[n4*64 + 42], arr_index[n4*64 + 46], n4%2==0);
        compareSwap(arr[n4*64 + 43], arr[n4*64 + 47], arr_index[n4*64 + 43], arr_index[n4*64 + 47], n4%2==0);

        compareSwap(arr[n4*64 + 48], arr[n4*64 + 52], arr_index[n4*64 + 48], arr_index[n4*64 + 52], n4%2==0);
        compareSwap(arr[n4*64 + 49], arr[n4*64 + 53], arr_index[n4*64 + 49], arr_index[n4*64 + 53], n4%2==0);
        compareSwap(arr[n4*64 + 50], arr[n4*64 + 54], arr_index[n4*64 + 50], arr_index[n4*64 + 54], n4%2==0);
        compareSwap(arr[n4*64 + 51], arr[n4*64 + 55], arr_index[n4*64 + 51], arr_index[n4*64 + 55], n4%2==0);

        compareSwap(arr[n4*64 + 56], arr[n4*64 + 60], arr_index[n4*64 + 56], arr_index[n4*64 + 60], n4%2==0);
        compareSwap(arr[n4*64 + 57], arr[n4*64 + 61], arr_index[n4*64 + 57], arr_index[n4*64 + 61], n4%2==0);
        compareSwap(arr[n4*64 + 58], arr[n4*64 + 62], arr_index[n4*64 + 58], arr_index[n4*64 + 62], n4%2==0);
        compareSwap(arr[n4*64 + 59], arr[n4*64 + 63], arr_index[n4*64 + 59], arr_index[n4*64 + 63], n4%2==0);
    }

    for(int n5 = 0; n5 < 4; n5++){

        compareSwap(arr[n5*64], arr[n5*64 + 2], arr_index[n5*64], arr_index[n5*64 + 2], n5%2==0);
        compareSwap(arr[n5*64 + 1], arr[n5*64 + 3], arr_index[n5*64 + 1], arr_index[n5*64 + 3], n5%2==0);
        compareSwap(arr[n5*64 + 4], arr[n5*64 + 6], arr_index[n5*64 + 4], arr_index[n5*64 + 6], n5%2==0);
        compareSwap(arr[n5*64 + 5], arr[n5*64 + 7], arr_index[n5*64 + 5], arr_index[n5*64 + 7], n5%2==0);
        compareSwap(arr[n5*64 + 8], arr[n5*64 + 10], arr_index[n5*64 + 8], arr_index[n5*64 + 10], n5%2==0);
        compareSwap(arr[n5*64 + 9], arr[n5*64 + 11], arr_index[n5*64 + 9], arr_index[n5*64 + 11], n5%2==0);
        compareSwap(arr[n5*64 + 12], arr[n5*64 + 14], arr_index[n5*64 + 12], arr_index[n5*64 + 14], n5%2==0);
        compareSwap(arr[n5*64 + 13], arr[n5*64 + 15], arr_index[n5*64 + 13], arr_index[n5*64 + 15], n5%2==0);

        compareSwap(arr[n5*64 + 16], arr[n5*64 + 18], arr_index[n5*64 + 16], arr_index[n5*64 + 18], n5%2==0);
        compareSwap(arr[n5*64 + 17], arr[n5*64 + 19], arr_index[n5*64 + 17], arr_index[n5*64 + 19], n5%2==0);
        compareSwap(arr[n5*64 + 20], arr[n5*64 + 22], arr_index[n5*64 + 20], arr_index[n5*64 + 22], n5%2==0);
        compareSwap(arr[n5*64 + 21], arr[n5*64 + 23], arr_index[n5*64 + 21], arr_index[n5*64 + 23], n5%2==0);
        compareSwap(arr[n5*64 + 24], arr[n5*64 + 26], arr_index[n5*64 + 24], arr_index[n5*64 + 26], n5%2==0);
        compareSwap(arr[n5*64 + 25], arr[n5*64 + 27], arr_index[n5*64 + 25], arr_index[n5*64 + 27], n5%2==0);
        compareSwap(arr[n5*64 + 28], arr[n5*64 + 30], arr_index[n5*64 + 28], arr_index[n5*64 + 30], n5%2==0);
        compareSwap(arr[n5*64 + 29], arr[n5*64 + 31], arr_index[n5*64 + 29], arr_index[n5*64 + 31], n5%2==0);
        compareSwap(arr[n5*64 + 32], arr[n5*64 + 34], arr_index[n5*64 + 32], arr_index[n5*64 + 34], n5%2==0);
        compareSwap(arr[n5*64 + 33], arr[n5*64 + 35], arr_index[n5*64 + 33], arr_index[n5*64 + 35], n5%2==0);
        
        compareSwap(arr[n5*64 + 36], arr[n5*64 + 38], arr_index[n5*64 + 36], arr_index[n5*64 + 38], n5%2==0);
        compareSwap(arr[n5*64 + 37], arr[n5*64 + 39], arr_index[n5*64 + 37], arr_index[n5*64 + 39], n5%2==0);
        compareSwap(arr[n5*64 + 40], arr[n5*64 + 42], arr_index[n5*64 + 40], arr_index[n5*64 + 42], n5%2==0);
        compareSwap(arr[n5*64 + 41], arr[n5*64 + 43], arr_index[n5*64 + 41], arr_index[n5*64 + 43], n5%2==0);
        compareSwap(arr[n5*64 + 44], arr[n5*64 + 46], arr_index[n5*64 + 44], arr_index[n5*64 + 46], n5%2==0);
        compareSwap(arr[n5*64 + 45], arr[n5*64 + 47], arr_index[n5*64 + 45], arr_index[n5*64 + 47], n5%2==0);
        compareSwap(arr[n5*64 + 48], arr[n5*64 + 50], arr_index[n5*64 + 48], arr_index[n5*64 + 50], n5%2==0);
        compareSwap(arr[n5*64 + 49], arr[n5*64 + 51], arr_index[n5*64 + 49], arr_index[n5*64 + 51], n5%2==0);
        compareSwap(arr[n5*64 + 52], arr[n5*64 + 54], arr_index[n5*64 + 52], arr_index[n5*64 + 54], n5%2==0);
        compareSwap(arr[n5*64 + 53], arr[n5*64 + 55], arr_index[n5*64 + 53], arr_index[n5*64 + 55], n5%2==0);
        compareSwap(arr[n5*64 + 56], arr[n5*64 + 58], arr_index[n5*64 + 56], arr_index[n5*64 + 58], n5%2==0);
        compareSwap(arr[n5*64 + 57], arr[n5*64 + 59], arr_index[n5*64 + 57], arr_index[n5*64 + 59], n5%2==0);
        compareSwap(arr[n5*64 + 60], arr[n5*64 + 62], arr_index[n5*64 + 60], arr_index[n5*64 + 62], n5%2==0);
        compareSwap(arr[n5*64 + 61], arr[n5*64 + 63], arr_index[n5*64 + 61], arr_index[n5*64 + 63], n5%2==0);
    }

    for(int n6 = 0; n6 < 4; n6++){

        compareSwap(arr[n6*64], arr[n6*64 + 1], arr_index[n6*64], arr_index[n6*64 + 1], n6%2==0);
        compareSwap(arr[n6*64 + 2], arr[n6*64 + 3], arr_index[n6*64 + 2], arr_index[n6*64 + 3], n6%2==0);
        compareSwap(arr[n6*64 + 4], arr[n6*64 + 5], arr_index[n6*64 + 4], arr_index[n6*64 + 5], n6%2==0);
        compareSwap(arr[n6*64 + 6], arr[n6*64 + 7], arr_index[n6*64 + 6], arr_index[n6*64 + 7], n6%2==0);
        compareSwap(arr[n6*64 + 8], arr[n6*64 + 9], arr_index[n6*64 + 8], arr_index[n6*64 + 9], n6%2==0);
        compareSwap(arr[n6*64 + 10], arr[n6*64 + 11], arr_index[n6*64 + 10], arr_index[n6*64 + 11], n6%2==0);
        compareSwap(arr[n6*64 + 12], arr[n6*64 + 13], arr_index[n6*64 + 12], arr_index[n6*64 + 13], n6%2==0);
        compareSwap(arr[n6*64 + 14], arr[n6*64 + 15], arr_index[n6*64 + 14], arr_index[n6*64 + 15], n6%2==0);
        compareSwap(arr[n6*64 + 16], arr[n6*64 + 17], arr_index[n6*64 + 16], arr_index[n6*64 + 17], n6%2==0);
        compareSwap(arr[n6*64 + 18], arr[n6*64 + 19], arr_index[n6*64 + 18], arr_index[n6*64 + 19], n6%2==0);
        compareSwap(arr[n6*64 + 20], arr[n6*64 + 21], arr_index[n6*64 + 20], arr_index[n6*64 + 21], n6%2==0);
        compareSwap(arr[n6*64 + 22], arr[n6*64 + 23], arr_index[n6*64 + 22], arr_index[n6*64 + 23], n6%2==0);
        compareSwap(arr[n6*64 + 24], arr[n6*64 + 25], arr_index[n6*64 + 24], arr_index[n6*64 + 25], n6%2==0);
        compareSwap(arr[n6*64 + 26], arr[n6*64 + 27], arr_index[n6*64 + 26], arr_index[n6*64 + 27], n6%2==0);
        compareSwap(arr[n6*64 + 28], arr[n6*64 + 29], arr_index[n6*64 + 28], arr_index[n6*64 + 29], n6%2==0);
        compareSwap(arr[n6*64 + 30], arr[n6*64 + 31], arr_index[n6*64 + 30], arr_index[n6*64 + 31], n6%2==0);
        compareSwap(arr[n6*64 + 32], arr[n6*64 + 33], arr_index[n6*64 + 32], arr_index[n6*64 + 33], n6%2==0);
        compareSwap(arr[n6*64 + 34], arr[n6*64 + 35], arr_index[n6*64 + 34], arr_index[n6*64 + 35], n6%2==0);
        compareSwap(arr[n6*64 + 36], arr[n6*64 + 37], arr_index[n6*64 + 36], arr_index[n6*64 + 37], n6%2==0);
        compareSwap(arr[n6*64 + 38], arr[n6*64 + 39], arr_index[n6*64 + 38], arr_index[n6*64 + 39], n6%2==0);
        compareSwap(arr[n6*64 + 40], arr[n6*64 + 41], arr_index[n6*64 + 40], arr_index[n6*64 + 41], n6%2==0);
        compareSwap(arr[n6*64 + 42], arr[n6*64 + 43], arr_index[n6*64 + 42], arr_index[n6*64 + 43], n6%2==0);
        compareSwap(arr[n6*64 + 44], arr[n6*64 + 45], arr_index[n6*64 + 44], arr_index[n6*64 + 45], n6%2==0);
        compareSwap(arr[n6*64 + 46], arr[n6*64 + 47], arr_index[n6*64 + 46], arr_index[n6*64 + 47], n6%2==0);
        compareSwap(arr[n6*64 + 48], arr[n6*64 + 49], arr_index[n6*64 + 48], arr_index[n6*64 + 49], n6%2==0);
        compareSwap(arr[n6*64 + 50], arr[n6*64 + 51], arr_index[n6*64 + 50], arr_index[n6*64 + 51], n6%2==0);
        compareSwap(arr[n6*64 + 52], arr[n6*64 + 53], arr_index[n6*64 + 52], arr_index[n6*64 + 53], n6%2==0);
        compareSwap(arr[n6*64 + 54], arr[n6*64 + 55], arr_index[n6*64 + 54], arr_index[n6*64 + 55], n6%2==0);
        compareSwap(arr[n6*64 + 56], arr[n6*64 + 57], arr_index[n6*64 + 56], arr_index[n6*64 + 57], n6%2==0);
        compareSwap(arr[n6*64 + 58], arr[n6*64 + 59], arr_index[n6*64 + 58], arr_index[n6*64 + 59], n6%2==0);
        compareSwap(arr[n6*64 + 60], arr[n6*64 + 61], arr_index[n6*64 + 60], arr_index[n6*64 + 61], n6%2==0);
        compareSwap(arr[n6*64 + 62], arr[n6*64 + 63], arr_index[n6*64 + 62], arr_index[n6*64 + 63], n6%2==0);
    }
    //Stage 7
    for(int o = 0; o < 2; o++){

        compareSwap(arr[o*128], arr[o*128 + 64], arr_index[o*128], arr_index[o*128 + 64], o%2==0);
        compareSwap(arr[o*128 + 1], arr[o*128 + 65], arr_index[o*128 + 1], arr_index[o*128 + 65], o%2==0);
        compareSwap(arr[o*128 + 2], arr[o*128 + 66], arr_index[o*128 + 2], arr_index[o*128 + 66], o%2==0);
        compareSwap(arr[o*128 + 3], arr[o*128 + 67], arr_index[o*128 + 3], arr_index[o*128 + 67], o%2==0);
        compareSwap(arr[o*128 + 4], arr[o*128 + 68], arr_index[o*128 + 4], arr_index[o*128 + 68], o%2==0);
        compareSwap(arr[o*128 + 5], arr[o*128 + 69], arr_index[o*128 + 5], arr_index[o*128 + 69], o%2==0);
        compareSwap(arr[o*128 + 6], arr[o*128 + 70], arr_index[o*128 + 6], arr_index[o*128 + 70], o%2==0);
        compareSwap(arr[o*128 + 7], arr[o*128 + 71], arr_index[o*128 + 7], arr_index[o*128 + 71], o%2==0);
        compareSwap(arr[o*128 + 8], arr[o*128 + 72], arr_index[o*128 + 8], arr_index[o*128 + 72], o%2==0);
        compareSwap(arr[o*128 + 9], arr[o*128 + 73], arr_index[o*128 + 9], arr_index[o*128 + 73], o%2==0);
        compareSwap(arr[o*128 + 10], arr[o*128 + 74], arr_index[o*128 + 10], arr_index[o*128 + 74], o%2==0);
        compareSwap(arr[o*128 + 11], arr[o*128 + 75], arr_index[o*128 + 11], arr_index[o*128 + 75], o%2==0);
        compareSwap(arr[o*128 + 12], arr[o*128 + 76], arr_index[o*128 + 12], arr_index[o*128 + 76], o%2==0);
        compareSwap(arr[o*128 + 13], arr[o*128 + 77], arr_index[o*128 + 13], arr_index[o*128 + 77], o%2==0);
        compareSwap(arr[o*128 + 14], arr[o*128 + 78], arr_index[o*128 + 14], arr_index[o*128 + 78], o%2==0);
        compareSwap(arr[o*128 + 15], arr[o*128 + 79], arr_index[o*128 + 15], arr_index[o*128 + 79], o%2==0);
        compareSwap(arr[o*128 + 16], arr[o*128 + 80], arr_index[o*128 + 16], arr_index[o*128 + 80], o%2==0);
        compareSwap(arr[o*128 + 17], arr[o*128 + 81], arr_index[o*128 + 17], arr_index[o*128 + 81], o%2==0);
        compareSwap(arr[o*128 + 18], arr[o*128 + 82], arr_index[o*128 + 18], arr_index[o*128 + 82], o%2==0);
        compareSwap(arr[o*128 + 19], arr[o*128 + 83], arr_index[o*128 + 19], arr_index[o*128 + 83], o%2==0);
        compareSwap(arr[o*128 + 20], arr[o*128 + 84], arr_index[o*128 + 20], arr_index[o*128 + 84], o%2==0);
        compareSwap(arr[o*128 + 21], arr[o*128 + 85], arr_index[o*128 + 21], arr_index[o*128 + 85], o%2==0);
        compareSwap(arr[o*128 + 22], arr[o*128 + 86], arr_index[o*128 + 22], arr_index[o*128 + 86], o%2==0);
        compareSwap(arr[o*128 + 23], arr[o*128 + 87], arr_index[o*128 + 23], arr_index[o*128 + 87], o%2==0);
        compareSwap(arr[o*128 + 24], arr[o*128 + 88], arr_index[o*128 + 24], arr_index[o*128 + 88], o%2==0);
        compareSwap(arr[o*128 + 25], arr[o*128 + 89], arr_index[o*128 + 25], arr_index[o*128 + 89], o%2==0);
        compareSwap(arr[o*128 + 26], arr[o*128 + 90], arr_index[o*128 + 26], arr_index[o*128 + 90], o%2==0);
        compareSwap(arr[o*128 + 27], arr[o*128 + 91], arr_index[o*128 + 27], arr_index[o*128 + 91], o%2==0);
        compareSwap(arr[o*128 + 28], arr[o*128 + 92], arr_index[o*128 + 28], arr_index[o*128 + 92], o%2==0);
        compareSwap(arr[o*128 + 29], arr[o*128 + 93], arr_index[o*128 + 29], arr_index[o*128 + 93], o%2==0);
        compareSwap(arr[o*128 + 30], arr[o*128 + 94], arr_index[o*128 + 30], arr_index[o*128 + 94], o%2==0);
        compareSwap(arr[o*128 + 31], arr[o*128 + 95], arr_index[o*128 + 31], arr_index[o*128 + 95], o%2==0);
        compareSwap(arr[o*128 + 32], arr[o*128 + 96], arr_index[o*128 + 32], arr_index[o*128 + 96], o%2==0);
        compareSwap(arr[o*128 + 33], arr[o*128 + 97], arr_index[o*128 + 33], arr_index[o*128 + 97], o%2==0);
        compareSwap(arr[o*128 + 34], arr[o*128 + 98], arr_index[o*128 + 34], arr_index[o*128 + 98], o%2==0);
        compareSwap(arr[o*128 + 35], arr[o*128 + 99], arr_index[o*128 + 35], arr_index[o*128 + 99], o%2==0);
        compareSwap(arr[o*128 + 36], arr[o*128 + 100], arr_index[o*128 + 36], arr_index[o*128 + 100], o%2==0);
        compareSwap(arr[o*128 + 37], arr[o*128 + 101], arr_index[o*128 + 37], arr_index[o*128 + 101], o%2==0);
        compareSwap(arr[o*128 + 38], arr[o*128 + 102], arr_index[o*128 + 38], arr_index[o*128 + 102], o%2==0);
        compareSwap(arr[o*128 + 39], arr[o*128 + 103], arr_index[o*128 + 39], arr_index[o*128 + 103], o%2==0);
        compareSwap(arr[o*128 + 40], arr[o*128 + 104], arr_index[o*128 + 40], arr_index[o*128 + 104], o%2==0);
        compareSwap(arr[o*128 + 41], arr[o*128 + 105], arr_index[o*128 + 41], arr_index[o*128 + 105], o%2==0);
        compareSwap(arr[o*128 + 42], arr[o*128 + 106], arr_index[o*128 + 42], arr_index[o*128 + 106], o%2==0);
        compareSwap(arr[o*128 + 43], arr[o*128 + 107], arr_index[o*128 + 43], arr_index[o*128 + 107], o%2==0);
        compareSwap(arr[o*128 + 44], arr[o*128 + 108], arr_index[o*128 + 44], arr_index[o*128 + 108], o%2==0);
        compareSwap(arr[o*128 + 45], arr[o*128 + 109], arr_index[o*128 + 45], arr_index[o*128 + 109], o%2==0);
        compareSwap(arr[o*128 + 46], arr[o*128 + 110], arr_index[o*128 + 46], arr_index[o*128 + 110], o%2==0);
        compareSwap(arr[o*128 + 47], arr[o*128 + 111], arr_index[o*128 + 47], arr_index[o*128 + 111], o%2==0);
        compareSwap(arr[o*128 + 48], arr[o*128 + 112], arr_index[o*128 + 48], arr_index[o*128 + 112], o%2==0);
        compareSwap(arr[o*128 + 49], arr[o*128 + 113], arr_index[o*128 + 49], arr_index[o*128 + 113], o%2==0);
        compareSwap(arr[o*128 + 50], arr[o*128 + 114], arr_index[o*128 + 50], arr_index[o*128 + 114], o%2==0);
        compareSwap(arr[o*128 + 51], arr[o*128 + 115], arr_index[o*128 + 51], arr_index[o*128 + 115], o%2==0);
        compareSwap(arr[o*128 + 52], arr[o*128 + 116], arr_index[o*128 + 52], arr_index[o*128 + 116], o%2==0);
        compareSwap(arr[o*128 + 53], arr[o*128 + 117], arr_index[o*128 + 53], arr_index[o*128 + 117], o%2==0);
        compareSwap(arr[o*128 + 54], arr[o*128 + 118], arr_index[o*128 + 54], arr_index[o*128 + 118], o%2==0);
        compareSwap(arr[o*128 + 55], arr[o*128 + 119], arr_index[o*128 + 55], arr_index[o*128 + 119], o%2==0);
        compareSwap(arr[o*128 + 56], arr[o*128 + 120], arr_index[o*128 + 56], arr_index[o*128 + 120], o%2==0);
        compareSwap(arr[o*128 + 57], arr[o*128 + 121], arr_index[o*128 + 57], arr_index[o*128 + 121], o%2==0);
        compareSwap(arr[o*128 + 58], arr[o*128 + 122], arr_index[o*128 + 58], arr_index[o*128 + 122], o%2==0);
        compareSwap(arr[o*128 + 59], arr[o*128 + 123], arr_index[o*128 + 59], arr_index[o*128 + 123], o%2==0);
        compareSwap(arr[o*128 + 60], arr[o*128 + 124], arr_index[o*128 + 60], arr_index[o*128 + 124], o%2==0);
        compareSwap(arr[o*128 + 61], arr[o*128 + 125], arr_index[o*128 + 61], arr_index[o*128 + 125], o%2==0);
        compareSwap(arr[o*128 + 62], arr[o*128 + 126], arr_index[o*128 + 62], arr_index[o*128 + 126], o%2==0);
        compareSwap(arr[o*128 + 63], arr[o*128 + 127], arr_index[o*128 + 63], arr_index[o*128 + 127], o%2==0);
    }
    for(int o1 = 0; o1 < 2; o1++){

        compareSwap(arr[o1*128], arr[o1*128 + 32], arr_index[o1*128], arr_index[o1*128 + 32], o1%2==0);
        compareSwap(arr[o1*128 + 1], arr[o1*128 + 33], arr_index[o1*128 + 1], arr_index[o1*128 + 33], o1%2==0);
        compareSwap(arr[o1*128 + 2], arr[o1*128 + 34], arr_index[o1*128 + 2], arr_index[o1*128 + 34], o1%2==0);
        compareSwap(arr[o1*128 + 3], arr[o1*128 + 35], arr_index[o1*128 + 3], arr_index[o1*128 + 35], o1%2==0);
        compareSwap(arr[o1*128 + 4], arr[o1*128 + 36], arr_index[o1*128 + 4], arr_index[o1*128 + 36], o1%2==0);
        compareSwap(arr[o1*128 + 5], arr[o1*128 + 37], arr_index[o1*128 + 5], arr_index[o1*128 + 37], o1%2==0);
        compareSwap(arr[o1*128 + 6], arr[o1*128 + 38], arr_index[o1*128 + 6], arr_index[o1*128 + 38], o1%2==0);
        compareSwap(arr[o1*128 + 7], arr[o1*128 + 39], arr_index[o1*128 + 7], arr_index[o1*128 + 39], o1%2==0);
        compareSwap(arr[o1*128 + 8], arr[o1*128 + 40], arr_index[o1*128 + 8], arr_index[o1*128 + 40], o1%2==0);
        compareSwap(arr[o1*128 + 9], arr[o1*128 + 41], arr_index[o1*128 + 9], arr_index[o1*128 + 41], o1%2==0);
        compareSwap(arr[o1*128 + 10], arr[o1*128 + 42], arr_index[o1*128 + 10], arr_index[o1*128 + 42], o1%2==0);
        compareSwap(arr[o1*128 + 11], arr[o1*128 + 43], arr_index[o1*128 + 11], arr_index[o1*128 + 43], o1%2==0);
        compareSwap(arr[o1*128 + 12], arr[o1*128 + 44], arr_index[o1*128 + 12], arr_index[o1*128 + 44], o1%2==0);
        compareSwap(arr[o1*128 + 13], arr[o1*128 + 45], arr_index[o1*128 + 13], arr_index[o1*128 + 45], o1%2==0);
        compareSwap(arr[o1*128 + 14], arr[o1*128 + 46], arr_index[o1*128 + 14], arr_index[o1*128 + 46], o1%2==0);
        compareSwap(arr[o1*128 + 15], arr[o1*128 + 47], arr_index[o1*128 + 15], arr_index[o1*128 + 47], o1%2==0);
        compareSwap(arr[o1*128 + 16], arr[o1*128 + 48], arr_index[o1*128 + 16], arr_index[o1*128 + 48], o1%2==0);
        compareSwap(arr[o1*128 + 17], arr[o1*128 + 49], arr_index[o1*128 + 17], arr_index[o1*128 + 49], o1%2==0);
        compareSwap(arr[o1*128 + 18], arr[o1*128 + 50], arr_index[o1*128 + 18], arr_index[o1*128 + 50], o1%2==0);
        compareSwap(arr[o1*128 + 19], arr[o1*128 + 51], arr_index[o1*128 + 19], arr_index[o1*128 + 51], o1%2==0);
        compareSwap(arr[o1*128 + 20], arr[o1*128 + 52], arr_index[o1*128 + 20], arr_index[o1*128 + 52], o1%2==0);
        compareSwap(arr[o1*128 + 21], arr[o1*128 + 53], arr_index[o1*128 + 21], arr_index[o1*128 + 53], o1%2==0);
        compareSwap(arr[o1*128 + 22], arr[o1*128 + 54], arr_index[o1*128 + 22], arr_index[o1*128 + 54], o1%2==0);
        compareSwap(arr[o1*128 + 23], arr[o1*128 + 55], arr_index[o1*128 + 23], arr_index[o1*128 + 55], o1%2==0);
        compareSwap(arr[o1*128 + 24], arr[o1*128 + 56], arr_index[o1*128 + 24], arr_index[o1*128 + 56], o1%2==0);
        compareSwap(arr[o1*128 + 25], arr[o1*128 + 57], arr_index[o1*128 + 25], arr_index[o1*128 + 57], o1%2==0);
        compareSwap(arr[o1*128 + 26], arr[o1*128 + 58], arr_index[o1*128 + 26], arr_index[o1*128 + 58], o1%2==0);
        compareSwap(arr[o1*128 + 27], arr[o1*128 + 59], arr_index[o1*128 + 27], arr_index[o1*128 + 59], o1%2==0);
        compareSwap(arr[o1*128 + 28], arr[o1*128 + 60], arr_index[o1*128 + 28], arr_index[o1*128 + 60], o1%2==0);
        compareSwap(arr[o1*128 + 29], arr[o1*128 + 61], arr_index[o1*128 + 29], arr_index[o1*128 + 61], o1%2==0);
        compareSwap(arr[o1*128 + 30], arr[o1*128 + 62], arr_index[o1*128 + 30], arr_index[o1*128 + 62], o1%2==0);
        compareSwap(arr[o1*128 + 31], arr[o1*128 + 63], arr_index[o1*128 + 31], arr_index[o1*128 + 63], o1%2==0);
        
        compareSwap(arr[o1*128 + 64], arr[o1*128 + 96], arr_index[o1*128 + 64], arr_index[o1*128 + 96], o1%2==0);
        compareSwap(arr[o1*128 + 65], arr[o1*128 + 97], arr_index[o1*128 + 65], arr_index[o1*128 + 97], o1%2==0);
        compareSwap(arr[o1*128 + 66], arr[o1*128 + 98], arr_index[o1*128 + 66], arr_index[o1*128 + 98], o1%2==0);
        compareSwap(arr[o1*128 + 67], arr[o1*128 + 99], arr_index[o1*128 + 67], arr_index[o1*128 + 99], o1%2==0);
        compareSwap(arr[o1*128 + 68], arr[o1*128 + 100], arr_index[o1*128 + 68], arr_index[o1*128 + 100], o1%2==0);
        compareSwap(arr[o1*128 + 69], arr[o1*128 + 101], arr_index[o1*128 + 69], arr_index[o1*128 + 101], o1%2==0);
        compareSwap(arr[o1*128 + 70], arr[o1*128 + 102], arr_index[o1*128 + 70], arr_index[o1*128 + 102], o1%2==0);
        compareSwap(arr[o1*128 + 71], arr[o1*128 + 103], arr_index[o1*128 + 71], arr_index[o1*128 + 103], o1%2==0);
        compareSwap(arr[o1*128 + 72], arr[o1*128 + 104], arr_index[o1*128 + 72], arr_index[o1*128 + 104], o1%2==0);
        compareSwap(arr[o1*128 + 73], arr[o1*128 + 105], arr_index[o1*128 + 73], arr_index[o1*128 + 105], o1%2==0);
        compareSwap(arr[o1*128 + 74], arr[o1*128 + 106], arr_index[o1*128 + 74], arr_index[o1*128 + 106], o1%2==0);
        compareSwap(arr[o1*128 + 75], arr[o1*128 + 107], arr_index[o1*128 + 75], arr_index[o1*128 + 107], o1%2==0);
        compareSwap(arr[o1*128 + 76], arr[o1*128 + 108], arr_index[o1*128 + 76], arr_index[o1*128 + 108], o1%2==0);
        compareSwap(arr[o1*128 + 77], arr[o1*128 + 109], arr_index[o1*128 + 77], arr_index[o1*128 + 109], o1%2==0);
        compareSwap(arr[o1*128 + 78], arr[o1*128 + 110], arr_index[o1*128 + 78], arr_index[o1*128 + 110], o1%2==0);
        compareSwap(arr[o1*128 + 79], arr[o1*128 + 111], arr_index[o1*128 + 79], arr_index[o1*128 + 111], o1%2==0);
        compareSwap(arr[o1*128 + 80], arr[o1*128 + 112], arr_index[o1*128 + 80], arr_index[o1*128 + 112], o1%2==0);
        compareSwap(arr[o1*128 + 81], arr[o1*128 + 113], arr_index[o1*128 + 81], arr_index[o1*128 + 113], o1%2==0);
        compareSwap(arr[o1*128 + 82], arr[o1*128 + 114], arr_index[o1*128 + 82], arr_index[o1*128 + 114], o1%2==0);
        compareSwap(arr[o1*128 + 83], arr[o1*128 + 115], arr_index[o1*128 + 83], arr_index[o1*128 + 115], o1%2==0);
        compareSwap(arr[o1*128 + 84], arr[o1*128 + 116], arr_index[o1*128 + 84], arr_index[o1*128 + 116], o1%2==0);
        compareSwap(arr[o1*128 + 85], arr[o1*128 + 117], arr_index[o1*128 + 85], arr_index[o1*128 + 117], o1%2==0);
        compareSwap(arr[o1*128 + 86], arr[o1*128 + 118], arr_index[o1*128 + 86], arr_index[o1*128 + 118], o1%2==0);
        compareSwap(arr[o1*128 + 87], arr[o1*128 + 119], arr_index[o1*128 + 87], arr_index[o1*128 + 119], o1%2==0);
        compareSwap(arr[o1*128 + 88], arr[o1*128 + 120], arr_index[o1*128 + 88], arr_index[o1*128 + 120], o1%2==0);
        compareSwap(arr[o1*128 + 89], arr[o1*128 + 121], arr_index[o1*128 + 89], arr_index[o1*128 + 121], o1%2==0);
        compareSwap(arr[o1*128 + 90], arr[o1*128 + 122], arr_index[o1*128 + 90], arr_index[o1*128 + 122], o1%2==0);
        compareSwap(arr[o1*128 + 91], arr[o1*128 + 123], arr_index[o1*128 + 91], arr_index[o1*128 + 123], o1%2==0);
        compareSwap(arr[o1*128 + 92], arr[o1*128 + 124], arr_index[o1*128 + 92], arr_index[o1*128 + 124], o1%2==0);
        compareSwap(arr[o1*128 + 93], arr[o1*128 + 125], arr_index[o1*128 + 93], arr_index[o1*128 + 125], o1%2==0);
        compareSwap(arr[o1*128 + 94], arr[o1*128 + 126], arr_index[o1*128 + 94], arr_index[o1*128 + 126], o1%2==0);
        compareSwap(arr[o1*128 + 95], arr[o1*128 + 127], arr_index[o1*128 + 95], arr_index[o1*128 + 127], o1%2==0);
    }
    for(int o2 = 0; o2 < 2; o2++){

        compareSwap(arr[o2*128], arr[o2*128 + 16], arr_index[o2*128], arr_index[o2*128 + 16], o2%2==0);
        compareSwap(arr[o2*128 + 1], arr[o2*128 + 17], arr_index[o2*128 + 1], arr_index[o2*128 + 17], o2%2==0);
        compareSwap(arr[o2*128 + 2], arr[o2*128 + 18], arr_index[o2*128 + 2], arr_index[o2*128 + 18], o2%2==0);
        compareSwap(arr[o2*128 + 3], arr[o2*128 + 19], arr_index[o2*128 + 3], arr_index[o2*128 + 19], o2%2==0);
        compareSwap(arr[o2*128 + 4], arr[o2*128 + 20], arr_index[o2*128 + 4], arr_index[o2*128 + 20], o2%2==0);
        compareSwap(arr[o2*128 + 5], arr[o2*128 + 21], arr_index[o2*128 + 5], arr_index[o2*128 + 21], o2%2==0);
        compareSwap(arr[o2*128 + 6], arr[o2*128 + 22], arr_index[o2*128 + 6], arr_index[o2*128 + 22], o2%2==0);
        compareSwap(arr[o2*128 + 7], arr[o2*128 + 23], arr_index[o2*128 + 7], arr_index[o2*128 + 23], o2%2==0);
        compareSwap(arr[o2*128 + 8], arr[o2*128 + 24], arr_index[o2*128 + 8], arr_index[o2*128 + 24], o2%2==0);
        compareSwap(arr[o2*128 + 9], arr[o2*128 + 25], arr_index[o2*128 + 9], arr_index[o2*128 + 25], o2%2==0);
        compareSwap(arr[o2*128 + 10], arr[o2*128 + 26], arr_index[o2*128 + 10], arr_index[o2*128 + 26], o2%2==0);
        compareSwap(arr[o2*128 + 11], arr[o2*128 + 27], arr_index[o2*128 + 11], arr_index[o2*128 + 27], o2%2==0);
        compareSwap(arr[o2*128 + 12], arr[o2*128 + 28], arr_index[o2*128 + 12], arr_index[o2*128 + 28], o2%2==0);
        compareSwap(arr[o2*128 + 13], arr[o2*128 + 29], arr_index[o2*128 + 13], arr_index[o2*128 + 29], o2%2==0);
        compareSwap(arr[o2*128 + 14], arr[o2*128 + 30], arr_index[o2*128 + 14], arr_index[o2*128 + 30], o2%2==0);
        compareSwap(arr[o2*128 + 15], arr[o2*128 + 31], arr_index[o2*128 + 15], arr_index[o2*128 + 31], o2%2==0);

        compareSwap(arr[o2*128 + 32], arr[o2*128 + 48], arr_index[o2*128 + 32], arr_index[o2*128 + 48], o2%2==0);
        compareSwap(arr[o2*128 + 33], arr[o2*128 + 49], arr_index[o2*128 + 33], arr_index[o2*128 + 49], o2%2==0);
        compareSwap(arr[o2*128 + 34], arr[o2*128 + 50], arr_index[o2*128 + 34], arr_index[o2*128 + 50], o2%2==0);
        compareSwap(arr[o2*128 + 35], arr[o2*128 + 51], arr_index[o2*128 + 35], arr_index[o2*128 + 51], o2%2==0);
        compareSwap(arr[o2*128 + 36], arr[o2*128 + 52], arr_index[o2*128 + 36], arr_index[o2*128 + 52], o2%2==0);
        compareSwap(arr[o2*128 + 37], arr[o2*128 + 53], arr_index[o2*128 + 37], arr_index[o2*128 + 53], o2%2==0);
        compareSwap(arr[o2*128 + 38], arr[o2*128 + 54], arr_index[o2*128 + 38], arr_index[o2*128 + 54], o2%2==0);
        compareSwap(arr[o2*128 + 39], arr[o2*128 + 55], arr_index[o2*128 + 39], arr_index[o2*128 + 55], o2%2==0);
        compareSwap(arr[o2*128 + 40], arr[o2*128 + 56], arr_index[o2*128 + 40], arr_index[o2*128 + 56], o2%2==0);
        compareSwap(arr[o2*128 + 41], arr[o2*128 + 57], arr_index[o2*128 + 41], arr_index[o2*128 + 57], o2%2==0);
        compareSwap(arr[o2*128 + 42], arr[o2*128 + 58], arr_index[o2*128 + 42], arr_index[o2*128 + 58], o2%2==0);
        compareSwap(arr[o2*128 + 43], arr[o2*128 + 59], arr_index[o2*128 + 43], arr_index[o2*128 + 59], o2%2==0);
        compareSwap(arr[o2*128 + 44], arr[o2*128 + 60], arr_index[o2*128 + 44], arr_index[o2*128 + 60], o2%2==0);
        compareSwap(arr[o2*128 + 45], arr[o2*128 + 61], arr_index[o2*128 + 45], arr_index[o2*128 + 61], o2%2==0);
        compareSwap(arr[o2*128 + 46], arr[o2*128 + 62], arr_index[o2*128 + 46], arr_index[o2*128 + 62], o2%2==0);
        compareSwap(arr[o2*128 + 47], arr[o2*128 + 63], arr_index[o2*128 + 47], arr_index[o2*128 + 63], o2%2==0);

        compareSwap(arr[o2*128 + 64], arr[o2*128 + 80], arr_index[o2*128 + 64], arr_index[o2*128 + 80], o2%2==0);
        compareSwap(arr[o2*128 + 65], arr[o2*128 + 81], arr_index[o2*128 + 65], arr_index[o2*128 + 81], o2%2==0);
        compareSwap(arr[o2*128 + 66], arr[o2*128 + 82], arr_index[o2*128 + 66], arr_index[o2*128 + 82], o2%2==0);
        compareSwap(arr[o2*128 + 67], arr[o2*128 + 83], arr_index[o2*128 + 67], arr_index[o2*128 + 83], o2%2==0);
        compareSwap(arr[o2*128 + 68], arr[o2*128 + 84], arr_index[o2*128 + 68], arr_index[o2*128 + 84], o2%2==0);
        compareSwap(arr[o2*128 + 69], arr[o2*128 + 85], arr_index[o2*128 + 69], arr_index[o2*128 + 85], o2%2==0);
        compareSwap(arr[o2*128 + 70], arr[o2*128 + 86], arr_index[o2*128 + 70], arr_index[o2*128 + 86], o2%2==0);
        compareSwap(arr[o2*128 + 71], arr[o2*128 + 87], arr_index[o2*128 + 71], arr_index[o2*128 + 87], o2%2==0);
        compareSwap(arr[o2*128 + 72], arr[o2*128 + 88], arr_index[o2*128 + 72], arr_index[o2*128 + 88], o2%2==0);
        compareSwap(arr[o2*128 + 73], arr[o2*128 + 89], arr_index[o2*128 + 73], arr_index[o2*128 + 89], o2%2==0);
        compareSwap(arr[o2*128 + 74], arr[o2*128 + 90], arr_index[o2*128 + 74], arr_index[o2*128 + 90], o2%2==0);
        compareSwap(arr[o2*128 + 75], arr[o2*128 + 91], arr_index[o2*128 + 75], arr_index[o2*128 + 91], o2%2==0);
        compareSwap(arr[o2*128 + 76], arr[o2*128 + 92], arr_index[o2*128 + 76], arr_index[o2*128 + 92], o2%2==0);
        compareSwap(arr[o2*128 + 77], arr[o2*128 + 93], arr_index[o2*128 + 77], arr_index[o2*128 + 93], o2%2==0);
        compareSwap(arr[o2*128 + 78], arr[o2*128 + 94], arr_index[o2*128 + 78], arr_index[o2*128 + 94], o2%2==0);
        compareSwap(arr[o2*128 + 79], arr[o2*128 + 95], arr_index[o2*128 + 79], arr_index[o2*128 + 95], o2%2==0);

        compareSwap(arr[o2*128 + 96], arr[o2*128 + 112], arr_index[o2*128 + 96], arr_index[o2*128 + 112], o2%2==0);
        compareSwap(arr[o2*128 + 97], arr[o2*128 + 113], arr_index[o2*128 + 97], arr_index[o2*128 + 113], o2%2==0);
        compareSwap(arr[o2*128 + 98], arr[o2*128 + 114], arr_index[o2*128 + 98], arr_index[o2*128 + 114], o2%2==0);
        compareSwap(arr[o2*128 + 99], arr[o2*128 + 115], arr_index[o2*128 + 99], arr_index[o2*128 + 115], o2%2==0);
        compareSwap(arr[o2*128 + 100], arr[o2*128 + 116], arr_index[o2*128 + 100], arr_index[o2*128 + 116], o2%2==0);
        compareSwap(arr[o2*128 + 101], arr[o2*128 + 117], arr_index[o2*128 + 101], arr_index[o2*128 + 117], o2%2==0);
        compareSwap(arr[o2*128 + 102], arr[o2*128 + 118], arr_index[o2*128 + 102], arr_index[o2*128 + 118], o2%2==0);
        compareSwap(arr[o2*128 + 103], arr[o2*128 + 119], arr_index[o2*128 + 103], arr_index[o2*128 + 119], o2%2==0);
        compareSwap(arr[o2*128 + 104], arr[o2*128 + 120], arr_index[o2*128 + 104], arr_index[o2*128 + 120], o2%2==0);
        compareSwap(arr[o2*128 + 105], arr[o2*128 + 121], arr_index[o2*128 + 105], arr_index[o2*128 + 121], o2%2==0);
        compareSwap(arr[o2*128 + 106], arr[o2*128 + 122], arr_index[o2*128 + 106], arr_index[o2*128 + 122], o2%2==0);
        compareSwap(arr[o2*128 + 107], arr[o2*128 + 123], arr_index[o2*128 + 107], arr_index[o2*128 + 123], o2%2==0);
        compareSwap(arr[o2*128 + 108], arr[o2*128 + 124], arr_index[o2*128 + 108], arr_index[o2*128 + 124], o2%2==0);
        compareSwap(arr[o2*128 + 109], arr[o2*128 + 125], arr_index[o2*128 + 109], arr_index[o2*128 + 125], o2%2==0);
        compareSwap(arr[o2*128 + 110], arr[o2*128 + 126], arr_index[o2*128 + 110], arr_index[o2*128 + 126], o2%2==0);
        compareSwap(arr[o2*128 + 111], arr[o2*128 + 127], arr_index[o2*128 + 111], arr_index[o2*128 + 127], o2%2==0);
    }
    for(int o3 = 0; o3 < 2; o3++){

        compareSwap(arr[o3*128], arr[o3*128 + 8], arr_index[o3*128], arr_index[o3*128 + 8], o3%2==0);
        compareSwap(arr[o3*128 + 1], arr[o3*128 + 9], arr_index[o3*128 + 1], arr_index[o3*128 + 9], o3%2==0);
        compareSwap(arr[o3*128 + 2], arr[o3*128 + 10], arr_index[o3*128 + 2], arr_index[o3*128 + 10], o3%2==0);
        compareSwap(arr[o3*128 + 3], arr[o3*128 + 11], arr_index[o3*128 + 3], arr_index[o3*128 + 11], o3%2==0);
        compareSwap(arr[o3*128 + 4], arr[o3*128 + 12], arr_index[o3*128 + 4], arr_index[o3*128 + 12], o3%2==0);
        compareSwap(arr[o3*128 + 5], arr[o3*128 + 13], arr_index[o3*128 + 5], arr_index[o3*128 + 13], o3%2==0);
        compareSwap(arr[o3*128 + 6], arr[o3*128 + 14], arr_index[o3*128 + 6], arr_index[o3*128 + 14], o3%2==0);
        compareSwap(arr[o3*128 + 7], arr[o3*128 + 15], arr_index[o3*128 + 7], arr_index[o3*128 + 15], o3%2==0);

        compareSwap(arr[o3*128 + 16], arr[o3*128 + 24], arr_index[o3*128 + 16], arr_index[o3*128 + 24], o3%2==0);
        compareSwap(arr[o3*128 + 17], arr[o3*128 + 25], arr_index[o3*128 + 17], arr_index[o3*128 + 25], o3%2==0);
        compareSwap(arr[o3*128 + 18], arr[o3*128 + 26], arr_index[o3*128 + 18], arr_index[o3*128 + 26], o3%2==0);
        compareSwap(arr[o3*128 + 19], arr[o3*128 + 27], arr_index[o3*128 + 19], arr_index[o3*128 + 27], o3%2==0);
        compareSwap(arr[o3*128 + 20], arr[o3*128 + 28], arr_index[o3*128 + 20], arr_index[o3*128 + 28], o3%2==0);
        compareSwap(arr[o3*128 + 21], arr[o3*128 + 29], arr_index[o3*128 + 21], arr_index[o3*128 + 29], o3%2==0);
        compareSwap(arr[o3*128 + 22], arr[o3*128 + 30], arr_index[o3*128 + 22], arr_index[o3*128 + 30], o3%2==0);
        compareSwap(arr[o3*128 + 23], arr[o3*128 + 31], arr_index[o3*128 + 23], arr_index[o3*128 + 31], o3%2==0);

        compareSwap(arr[o3*128 + 32], arr[o3*128 + 40], arr_index[o3*128 + 32], arr_index[o3*128 + 40], o3%2==0);
        compareSwap(arr[o3*128 + 33], arr[o3*128 + 41], arr_index[o3*128 + 33], arr_index[o3*128 + 41], o3%2==0);
        compareSwap(arr[o3*128 + 34], arr[o3*128 + 42], arr_index[o3*128 + 34], arr_index[o3*128 + 42], o3%2==0);
        compareSwap(arr[o3*128 + 35], arr[o3*128 + 43], arr_index[o3*128 + 35], arr_index[o3*128 + 43], o3%2==0);
        compareSwap(arr[o3*128 + 36], arr[o3*128 + 44], arr_index[o3*128 + 36], arr_index[o3*128 + 44], o3%2==0);
        compareSwap(arr[o3*128 + 37], arr[o3*128 + 45], arr_index[o3*128 + 37], arr_index[o3*128 + 45], o3%2==0);
        compareSwap(arr[o3*128 + 38], arr[o3*128 + 46], arr_index[o3*128 + 38], arr_index[o3*128 + 46], o3%2==0);
        compareSwap(arr[o3*128 + 39], arr[o3*128 + 47], arr_index[o3*128 + 39], arr_index[o3*128 + 47], o3%2==0);

        compareSwap(arr[o3*128 + 48], arr[o3*128 + 56], arr_index[o3*128 + 48], arr_index[o3*128 + 56], o3%2==0);
        compareSwap(arr[o3*128 + 49], arr[o3*128 + 57], arr_index[o3*128 + 49], arr_index[o3*128 + 57], o3%2==0);
        compareSwap(arr[o3*128 + 50], arr[o3*128 + 58], arr_index[o3*128 + 50], arr_index[o3*128 + 58], o3%2==0);
        compareSwap(arr[o3*128 + 51], arr[o3*128 + 59], arr_index[o3*128 + 51], arr_index[o3*128 + 59], o3%2==0);
        compareSwap(arr[o3*128 + 52], arr[o3*128 + 60], arr_index[o3*128 + 52], arr_index[o3*128 + 60], o3%2==0);
        compareSwap(arr[o3*128 + 53], arr[o3*128 + 61], arr_index[o3*128 + 53], arr_index[o3*128 + 61], o3%2==0);
        compareSwap(arr[o3*128 + 54], arr[o3*128 + 62], arr_index[o3*128 + 54], arr_index[o3*128 + 62], o3%2==0);
        compareSwap(arr[o3*128 + 55], arr[o3*128 + 63], arr_index[o3*128 + 55], arr_index[o3*128 + 63], o3%2==0);

        compareSwap(arr[o3*128 + 64], arr[o3*128 + 72], arr_index[o3*128 + 64], arr_index[o3*128 + 72], o3%2==0);
        compareSwap(arr[o3*128 + 65], arr[o3*128 + 73], arr_index[o3*128 + 65], arr_index[o3*128 + 73], o3%2==0);
        compareSwap(arr[o3*128 + 66], arr[o3*128 + 74], arr_index[o3*128 + 66], arr_index[o3*128 + 74], o3%2==0);
        compareSwap(arr[o3*128 + 67], arr[o3*128 + 75], arr_index[o3*128 + 67], arr_index[o3*128 + 75], o3%2==0);
        compareSwap(arr[o3*128 + 68], arr[o3*128 + 76], arr_index[o3*128 + 68], arr_index[o3*128 + 76], o3%2==0);
        compareSwap(arr[o3*128 + 69], arr[o3*128 + 77], arr_index[o3*128 + 69], arr_index[o3*128 + 77], o3%2==0);
        compareSwap(arr[o3*128 + 70], arr[o3*128 + 78], arr_index[o3*128 + 70], arr_index[o3*128 + 78], o3%2==0);
        compareSwap(arr[o3*128 + 71], arr[o3*128 + 79], arr_index[o3*128 + 71], arr_index[o3*128 + 79], o3%2==0);

        compareSwap(arr[o3*128 + 80], arr[o3*128 + 88], arr_index[o3*128 + 80], arr_index[o3*128 + 88], o3%2==0);
        compareSwap(arr[o3*128 + 81], arr[o3*128 + 89], arr_index[o3*128 + 81], arr_index[o3*128 + 89], o3%2==0);
        compareSwap(arr[o3*128 + 82], arr[o3*128 + 90], arr_index[o3*128 + 82], arr_index[o3*128 + 90], o3%2==0);
        compareSwap(arr[o3*128 + 83], arr[o3*128 + 91], arr_index[o3*128 + 83], arr_index[o3*128 + 91], o3%2==0);
        compareSwap(arr[o3*128 + 84], arr[o3*128 + 92], arr_index[o3*128 + 84], arr_index[o3*128 + 92], o3%2==0);
        compareSwap(arr[o3*128 + 85], arr[o3*128 + 93], arr_index[o3*128 + 85], arr_index[o3*128 + 93], o3%2==0);
        compareSwap(arr[o3*128 + 86], arr[o3*128 + 94], arr_index[o3*128 + 86], arr_index[o3*128 + 94], o3%2==0);
        compareSwap(arr[o3*128 + 87], arr[o3*128 + 95], arr_index[o3*128 + 87], arr_index[o3*128 + 95], o3%2==0);

        compareSwap(arr[o3*128 + 96], arr[o3*128 + 104], arr_index[o3*128 + 96], arr_index[o3*128 + 104], o3%2==0);
        compareSwap(arr[o3*128 + 97], arr[o3*128 + 105], arr_index[o3*128 + 97], arr_index[o3*128 + 105], o3%2==0);
        compareSwap(arr[o3*128 + 98], arr[o3*128 + 106], arr_index[o3*128 + 98], arr_index[o3*128 + 106], o3%2==0);
        compareSwap(arr[o3*128 + 99], arr[o3*128 + 107], arr_index[o3*128 + 99], arr_index[o3*128 + 107], o3%2==0);
        compareSwap(arr[o3*128 + 100], arr[o3*128 + 108], arr_index[o3*128 + 100], arr_index[o3*128 + 108], o3%2==0);
        compareSwap(arr[o3*128 + 101], arr[o3*128 + 109], arr_index[o3*128 + 101], arr_index[o3*128 + 109], o3%2==0);
        compareSwap(arr[o3*128 + 102], arr[o3*128 + 110], arr_index[o3*128 + 102], arr_index[o3*128 + 110], o3%2==0);
        compareSwap(arr[o3*128 + 103], arr[o3*128 + 111], arr_index[o3*128 + 103], arr_index[o3*128 + 111], o3%2==0);

        compareSwap(arr[o3*128 + 112], arr[o3*128 + 120], arr_index[o3*128 + 112], arr_index[o3*128 + 120], o3%2==0);
        compareSwap(arr[o3*128 + 113], arr[o3*128 + 121], arr_index[o3*128 + 113], arr_index[o3*128 + 121], o3%2==0);
        compareSwap(arr[o3*128 + 114], arr[o3*128 + 122], arr_index[o3*128 + 114], arr_index[o3*128 + 122], o3%2==0);
        compareSwap(arr[o3*128 + 115], arr[o3*128 + 123], arr_index[o3*128 + 115], arr_index[o3*128 + 123], o3%2==0);
        compareSwap(arr[o3*128 + 116], arr[o3*128 + 124], arr_index[o3*128 + 116], arr_index[o3*128 + 124], o3%2==0);
        compareSwap(arr[o3*128 + 117], arr[o3*128 + 125], arr_index[o3*128 + 117], arr_index[o3*128 + 125], o3%2==0);
        compareSwap(arr[o3*128 + 118], arr[o3*128 + 126], arr_index[o3*128 + 118], arr_index[o3*128 + 126], o3%2==0);
        compareSwap(arr[o3*128 + 119], arr[o3*128 + 127], arr_index[o3*128 + 119], arr_index[o3*128 + 127], o3%2==0);
    }
    for(int o4 = 0; o4 < 2; o4++){

        compareSwap(arr[o4*128], arr[o4*128 + 4], arr_index[o4*128], arr_index[o4*128 + 4], o4%2==0);
        compareSwap(arr[o4*128 + 1], arr[o4*128 + 5], arr_index[o4*128 + 1], arr_index[o4*128 + 5], o4%2==0);
        compareSwap(arr[o4*128 + 2], arr[o4*128 + 6], arr_index[o4*128 + 2], arr_index[o4*128 + 6], o4%2==0);
        compareSwap(arr[o4*128 + 3], arr[o4*128 + 7], arr_index[o4*128 + 3], arr_index[o4*128 + 7], o4%2==0);

        compareSwap(arr[o4*128 + 8], arr[o4*128 + 12], arr_index[o4*128 + 8], arr_index[o4*128 + 12], o4%2==0);
        compareSwap(arr[o4*128 + 9], arr[o4*128 + 13], arr_index[o4*128 + 9], arr_index[o4*128 + 13], o4%2==0);
        compareSwap(arr[o4*128 + 10], arr[o4*128 + 14], arr_index[o4*128 + 10], arr_index[o4*128 + 14], o4%2==0);
        compareSwap(arr[o4*128 + 11], arr[o4*128 + 15], arr_index[o4*128 + 11], arr_index[o4*128 + 15], o4%2==0);

        compareSwap(arr[o4*128 + 16], arr[o4*128 + 20], arr_index[o4*128 + 16], arr_index[o4*128 + 20], o4%2==0);
        compareSwap(arr[o4*128 + 17], arr[o4*128 + 21], arr_index[o4*128 + 17], arr_index[o4*128 + 21], o4%2==0);
        compareSwap(arr[o4*128 + 18], arr[o4*128 + 22], arr_index[o4*128 + 18], arr_index[o4*128 + 22], o4%2==0);
        compareSwap(arr[o4*128 + 19], arr[o4*128 + 23], arr_index[o4*128 + 19], arr_index[o4*128 + 23], o4%2==0);

        compareSwap(arr[o4*128 + 24], arr[o4*128 + 28], arr_index[o4*128 + 24], arr_index[o4*128 + 28], o4%2==0);
        compareSwap(arr[o4*128 + 25], arr[o4*128 + 29], arr_index[o4*128 + 25], arr_index[o4*128 + 29], o4%2==0);
        compareSwap(arr[o4*128 + 26], arr[o4*128 + 30], arr_index[o4*128 + 26], arr_index[o4*128 + 30], o4%2==0);
        compareSwap(arr[o4*128 + 27], arr[o4*128 + 31], arr_index[o4*128 + 27], arr_index[o4*128 + 31], o4%2==0);

        compareSwap(arr[o4*128 + 32], arr[o4*128 + 36], arr_index[o4*128 + 32], arr_index[o4*128 + 36], o4%2==0);
        compareSwap(arr[o4*128 + 33], arr[o4*128 + 37], arr_index[o4*128 + 33], arr_index[o4*128 + 37], o4%2==0);
        compareSwap(arr[o4*128 + 34], arr[o4*128 + 38], arr_index[o4*128 + 34], arr_index[o4*128 + 38], o4%2==0);
        compareSwap(arr[o4*128 + 35], arr[o4*128 + 39], arr_index[o4*128 + 35], arr_index[o4*128 + 39], o4%2==0);

        compareSwap(arr[o4*128 + 40], arr[o4*128 + 44], arr_index[o4*128 + 40], arr_index[o4*128 + 44], o4%2==0);
        compareSwap(arr[o4*128 + 41], arr[o4*128 + 45], arr_index[o4*128 + 41], arr_index[o4*128 + 45], o4%2==0);
        compareSwap(arr[o4*128 + 42], arr[o4*128 + 46], arr_index[o4*128 + 42], arr_index[o4*128 + 46], o4%2==0);
        compareSwap(arr[o4*128 + 43], arr[o4*128 + 47], arr_index[o4*128 + 43], arr_index[o4*128 + 47], o4%2==0);

        compareSwap(arr[o4*128 + 48], arr[o4*128 + 52], arr_index[o4*128 + 48], arr_index[o4*128 + 52], o4%2==0);
        compareSwap(arr[o4*128 + 49], arr[o4*128 + 53], arr_index[o4*128 + 49], arr_index[o4*128 + 53], o4%2==0);
        compareSwap(arr[o4*128 + 50], arr[o4*128 + 54], arr_index[o4*128 + 50], arr_index[o4*128 + 54], o4%2==0);
        compareSwap(arr[o4*128 + 51], arr[o4*128 + 55], arr_index[o4*128 + 51], arr_index[o4*128 + 55], o4%2==0);

        compareSwap(arr[o4*128 + 56], arr[o4*128 + 60], arr_index[o4*128 + 56], arr_index[o4*128 + 60], o4%2==0);
        compareSwap(arr[o4*128 + 57], arr[o4*128 + 61], arr_index[o4*128 + 57], arr_index[o4*128 + 61], o4%2==0);
        compareSwap(arr[o4*128 + 58], arr[o4*128 + 62], arr_index[o4*128 + 58], arr_index[o4*128 + 62], o4%2==0);
        compareSwap(arr[o4*128 + 59], arr[o4*128 + 63], arr_index[o4*128 + 59], arr_index[o4*128 + 63], o4%2==0);

        compareSwap(arr[o4*128 + 64], arr[o4*128 + 68], arr_index[o4*128 + 64], arr_index[o4*128 + 68], o4%2==0);
        compareSwap(arr[o4*128 + 65], arr[o4*128 + 69], arr_index[o4*128 + 65], arr_index[o4*128 + 69], o4%2==0);
        compareSwap(arr[o4*128 + 66], arr[o4*128 + 70], arr_index[o4*128 + 66], arr_index[o4*128 + 70], o4%2==0);
        compareSwap(arr[o4*128 + 67], arr[o4*128 + 71], arr_index[o4*128 + 67], arr_index[o4*128 + 71], o4%2==0);

        compareSwap(arr[o4*128 + 72], arr[o4*128 + 76], arr_index[o4*128 + 72], arr_index[o4*128 + 76], o4%2==0);
        compareSwap(arr[o4*128 + 73], arr[o4*128 + 77], arr_index[o4*128 + 73], arr_index[o4*128 + 77], o4%2==0);
        compareSwap(arr[o4*128 + 74], arr[o4*128 + 78], arr_index[o4*128 + 74], arr_index[o4*128 + 78], o4%2==0);
        compareSwap(arr[o4*128 + 75], arr[o4*128 + 79], arr_index[o4*128 + 75], arr_index[o4*128 + 79], o4%2==0);

        compareSwap(arr[o4*128 + 80], arr[o4*128 + 84], arr_index[o4*128 + 80], arr_index[o4*128 + 84], o4%2==0);
        compareSwap(arr[o4*128 + 81], arr[o4*128 + 85], arr_index[o4*128 + 81], arr_index[o4*128 + 85], o4%2==0);
        compareSwap(arr[o4*128 + 82], arr[o4*128 + 86], arr_index[o4*128 + 82], arr_index[o4*128 + 86], o4%2==0);
        compareSwap(arr[o4*128 + 83], arr[o4*128 + 87], arr_index[o4*128 + 83], arr_index[o4*128 + 87], o4%2==0);

        compareSwap(arr[o4*128 + 88], arr[o4*128 + 92], arr_index[o4*128 + 88], arr_index[o4*128 + 92], o4%2==0);
        compareSwap(arr[o4*128 + 89], arr[o4*128 + 93], arr_index[o4*128 + 89], arr_index[o4*128 + 93], o4%2==0);
        compareSwap(arr[o4*128 + 90], arr[o4*128 + 94], arr_index[o4*128 + 90], arr_index[o4*128 + 94], o4%2==0);
        compareSwap(arr[o4*128 + 91], arr[o4*128 + 95], arr_index[o4*128 + 91], arr_index[o4*128 + 95], o4%2==0);

        compareSwap(arr[o4*128 + 96], arr[o4*128 + 100], arr_index[o4*128 + 96], arr_index[o4*128 + 100], o4%2==0);
        compareSwap(arr[o4*128 + 97], arr[o4*128 + 101], arr_index[o4*128 + 97], arr_index[o4*128 + 101], o4%2==0);
        compareSwap(arr[o4*128 + 98], arr[o4*128 + 102], arr_index[o4*128 + 98], arr_index[o4*128 + 102], o4%2==0);
        compareSwap(arr[o4*128 + 99], arr[o4*128 + 103], arr_index[o4*128 + 99], arr_index[o4*128 + 103], o4%2==0);

        compareSwap(arr[o4*128 + 104], arr[o4*128 + 108], arr_index[o4*128 + 104], arr_index[o4*128 + 108], o4%2==0);
        compareSwap(arr[o4*128 + 105], arr[o4*128 + 109], arr_index[o4*128 + 105], arr_index[o4*128 + 109], o4%2==0);
        compareSwap(arr[o4*128 + 106], arr[o4*128 + 110], arr_index[o4*128 + 106], arr_index[o4*128 + 110], o4%2==0);
        compareSwap(arr[o4*128 + 107], arr[o4*128 + 111], arr_index[o4*128 + 107], arr_index[o4*128 + 111], o4%2==0);

        compareSwap(arr[o4*128 + 112], arr[o4*128 + 116], arr_index[o4*128 + 112], arr_index[o4*128 + 116], o4%2==0);
        compareSwap(arr[o4*128 + 113], arr[o4*128 + 117], arr_index[o4*128 + 113], arr_index[o4*128 + 117], o4%2==0);
        compareSwap(arr[o4*128 + 114], arr[o4*128 + 118], arr_index[o4*128 + 114], arr_index[o4*128 + 118], o4%2==0);
        compareSwap(arr[o4*128 + 115], arr[o4*128 + 119], arr_index[o4*128 + 115], arr_index[o4*128 + 119], o4%2==0);

        compareSwap(arr[o4*128 + 120], arr[o4*128 + 124], arr_index[o4*128 + 120], arr_index[o4*128 + 124], o4%2==0);
        compareSwap(arr[o4*128 + 121], arr[o4*128 + 125], arr_index[o4*128 + 121], arr_index[o4*128 + 125], o4%2==0);
        compareSwap(arr[o4*128 + 122], arr[o4*128 + 126], arr_index[o4*128 + 122], arr_index[o4*128 + 126], o4%2==0);
        compareSwap(arr[o4*128 + 123], arr[o4*128 + 127], arr_index[o4*128 + 123], arr_index[o4*128 + 127], o4%2==0);
    }
    for(int o5 = 0; o5 < 2; o5++){

        compareSwap(arr[o5*128], arr[o5*128 + 2], arr_index[o5*128], arr_index[o5*128 + 2], o5%2==0);
        compareSwap(arr[o5*128 + 1], arr[o5*128 + 3], arr_index[o5*128 + 1], arr_index[o5*128 + 3], o5%2==0);
        compareSwap(arr[o5*128 + 4], arr[o5*128 + 6], arr_index[o5*128 + 4], arr_index[o5*128 + 6], o5%2==0);
        compareSwap(arr[o5*128 + 5], arr[o5*128 + 7], arr_index[o5*128 + 5], arr_index[o5*128 + 7], o5%2==0);
        compareSwap(arr[o5*128 + 8], arr[o5*128 + 10], arr_index[o5*128 + 8], arr_index[o5*128 + 10], o5%2==0);
        compareSwap(arr[o5*128 + 9], arr[o5*128 + 11], arr_index[o5*128 + 9], arr_index[o5*128 + 11], o5%2==0);
        compareSwap(arr[o5*128 + 12], arr[o5*128 + 14], arr_index[o5*128 + 12], arr_index[o5*128 + 14], o5%2==0);
        compareSwap(arr[o5*128 + 13], arr[o5*128 + 15], arr_index[o5*128 + 13], arr_index[o5*128 + 15], o5%2==0);
        compareSwap(arr[o5*128 + 16], arr[o5*128 + 18], arr_index[o5*128 + 16], arr_index[o5*128 + 18], o5%2==0);
        compareSwap(arr[o5*128 + 17], arr[o5*128 + 19], arr_index[o5*128 + 17], arr_index[o5*128 + 19], o5%2==0);
        compareSwap(arr[o5*128 + 20], arr[o5*128 + 22], arr_index[o5*128 + 20], arr_index[o5*128 + 22], o5%2==0);
        compareSwap(arr[o5*128 + 21], arr[o5*128 + 23], arr_index[o5*128 + 21], arr_index[o5*128 + 23], o5%2==0);
        compareSwap(arr[o5*128 + 24], arr[o5*128 + 26], arr_index[o5*128 + 24], arr_index[o5*128 + 26], o5%2==0);
        compareSwap(arr[o5*128 + 25], arr[o5*128 + 27], arr_index[o5*128 + 25], arr_index[o5*128 + 27], o5%2==0);
        compareSwap(arr[o5*128 + 28], arr[o5*128 + 30], arr_index[o5*128 + 28], arr_index[o5*128 + 30], o5%2==0);
        compareSwap(arr[o5*128 + 29], arr[o5*128 + 31], arr_index[o5*128 + 29], arr_index[o5*128 + 31], o5%2==0);
        compareSwap(arr[o5*128 + 32], arr[o5*128 + 34], arr_index[o5*128 + 32], arr_index[o5*128 + 34], o5%2==0);
        compareSwap(arr[o5*128 + 33], arr[o5*128 + 35], arr_index[o5*128 + 33], arr_index[o5*128 + 35], o5%2==0);
        compareSwap(arr[o5*128 + 36], arr[o5*128 + 38], arr_index[o5*128 + 36], arr_index[o5*128 + 38], o5%2==0);
        compareSwap(arr[o5*128 + 37], arr[o5*128 + 39], arr_index[o5*128 + 37], arr_index[o5*128 + 39], o5%2==0);
        compareSwap(arr[o5*128 + 40], arr[o5*128 + 42], arr_index[o5*128 + 40], arr_index[o5*128 + 42], o5%2==0);
        compareSwap(arr[o5*128 + 41], arr[o5*128 + 43], arr_index[o5*128 + 41], arr_index[o5*128 + 43], o5%2==0);
        compareSwap(arr[o5*128 + 44], arr[o5*128 + 46], arr_index[o5*128 + 44], arr_index[o5*128 + 46], o5%2==0);
        compareSwap(arr[o5*128 + 45], arr[o5*128 + 47], arr_index[o5*128 + 45], arr_index[o5*128 + 47], o5%2==0);
        compareSwap(arr[o5*128 + 48], arr[o5*128 + 50], arr_index[o5*128 + 48], arr_index[o5*128 + 50], o5%2==0);
        compareSwap(arr[o5*128 + 49], arr[o5*128 + 51], arr_index[o5*128 + 49], arr_index[o5*128 + 51], o5%2==0);
        compareSwap(arr[o5*128 + 52], arr[o5*128 + 54], arr_index[o5*128 + 52], arr_index[o5*128 + 54], o5%2==0);
        compareSwap(arr[o5*128 + 53], arr[o5*128 + 55], arr_index[o5*128 + 53], arr_index[o5*128 + 55], o5%2==0);
        compareSwap(arr[o5*128 + 56], arr[o5*128 + 58], arr_index[o5*128 + 56], arr_index[o5*128 + 58], o5%2==0);
        compareSwap(arr[o5*128 + 57], arr[o5*128 + 59], arr_index[o5*128 + 57], arr_index[o5*128 + 59], o5%2==0);
        compareSwap(arr[o5*128 + 60], arr[o5*128 + 62], arr_index[o5*128 + 60], arr_index[o5*128 + 62], o5%2==0);
        compareSwap(arr[o5*128 + 61], arr[o5*128 + 63], arr_index[o5*128 + 61], arr_index[o5*128 + 63], o5%2==0);
        compareSwap(arr[o5*128 + 64], arr[o5*128 + 66], arr_index[o5*128 + 64], arr_index[o5*128 + 66], o5%2==0);
        compareSwap(arr[o5*128 + 65], arr[o5*128 + 67], arr_index[o5*128 + 65], arr_index[o5*128 + 67], o5%2==0);
        compareSwap(arr[o5*128 + 68], arr[o5*128 + 70], arr_index[o5*128 + 68], arr_index[o5*128 + 70], o5%2==0);
        compareSwap(arr[o5*128 + 69], arr[o5*128 + 71], arr_index[o5*128 + 69], arr_index[o5*128 + 71], o5%2==0);
        compareSwap(arr[o5*128 + 72], arr[o5*128 + 74], arr_index[o5*128 + 72], arr_index[o5*128 + 74], o5%2==0);
        compareSwap(arr[o5*128 + 73], arr[o5*128 + 75], arr_index[o5*128 + 73], arr_index[o5*128 + 75], o5%2==0);
        compareSwap(arr[o5*128 + 76], arr[o5*128 + 78], arr_index[o5*128 + 76], arr_index[o5*128 + 78], o5%2==0);
        compareSwap(arr[o5*128 + 77], arr[o5*128 + 79], arr_index[o5*128 + 77], arr_index[o5*128 + 79], o5%2==0);
        compareSwap(arr[o5*128 + 80], arr[o5*128 + 82], arr_index[o5*128 + 80], arr_index[o5*128 + 82], o5%2==0);
        compareSwap(arr[o5*128 + 81], arr[o5*128 + 83], arr_index[o5*128 + 81], arr_index[o5*128 + 83], o5%2==0);
        compareSwap(arr[o5*128 + 84], arr[o5*128 + 86], arr_index[o5*128 + 84], arr_index[o5*128 + 86], o5%2==0);
        compareSwap(arr[o5*128 + 85], arr[o5*128 + 87], arr_index[o5*128 + 85], arr_index[o5*128 + 87], o5%2==0);
        compareSwap(arr[o5*128 + 88], arr[o5*128 + 90], arr_index[o5*128 + 88], arr_index[o5*128 + 90], o5%2==0);
        compareSwap(arr[o5*128 + 89], arr[o5*128 + 91], arr_index[o5*128 + 89], arr_index[o5*128 + 91], o5%2==0);
        compareSwap(arr[o5*128 + 92], arr[o5*128 + 94], arr_index[o5*128 + 92], arr_index[o5*128 + 94], o5%2==0);
        compareSwap(arr[o5*128 + 93], arr[o5*128 + 95], arr_index[o5*128 + 93], arr_index[o5*128 + 95], o5%2==0);
        compareSwap(arr[o5*128 + 96], arr[o5*128 + 98], arr_index[o5*128 + 96], arr_index[o5*128 + 98], o5%2==0);
        compareSwap(arr[o5*128 + 97], arr[o5*128 + 99], arr_index[o5*128 + 97], arr_index[o5*128 + 99], o5%2==0);
        compareSwap(arr[o5*128 + 100], arr[o5*128 + 102], arr_index[o5*128 + 100], arr_index[o5*128 + 102], o5%2==0);
        compareSwap(arr[o5*128 + 101], arr[o5*128 + 103], arr_index[o5*128 + 101], arr_index[o5*128 + 103], o5%2==0);
        compareSwap(arr[o5*128 + 104], arr[o5*128 + 106], arr_index[o5*128 + 104], arr_index[o5*128 + 106], o5%2==0);
        compareSwap(arr[o5*128 + 105], arr[o5*128 + 107], arr_index[o5*128 + 105], arr_index[o5*128 + 107], o5%2==0);
        compareSwap(arr[o5*128 + 108], arr[o5*128 + 110], arr_index[o5*128 + 108], arr_index[o5*128 + 110], o5%2==0);
        compareSwap(arr[o5*128 + 109], arr[o5*128 + 111], arr_index[o5*128 + 109], arr_index[o5*128 + 111], o5%2==0);
        compareSwap(arr[o5*128 + 112], arr[o5*128 + 114], arr_index[o5*128 + 112], arr_index[o5*128 + 114], o5%2==0);
        compareSwap(arr[o5*128 + 113], arr[o5*128 + 115], arr_index[o5*128 + 113], arr_index[o5*128 + 115], o5%2==0);
        compareSwap(arr[o5*128 + 116], arr[o5*128 + 118], arr_index[o5*128 + 116], arr_index[o5*128 + 118], o5%2==0);
        compareSwap(arr[o5*128 + 117], arr[o5*128 + 119], arr_index[o5*128 + 117], arr_index[o5*128 + 119], o5%2==0);
        compareSwap(arr[o5*128 + 120], arr[o5*128 + 122], arr_index[o5*128 + 120], arr_index[o5*128 + 122], o5%2==0);
        compareSwap(arr[o5*128 + 121], arr[o5*128 + 123], arr_index[o5*128 + 121], arr_index[o5*128 + 123], o5%2==0);
        compareSwap(arr[o5*128 + 124], arr[o5*128 + 126], arr_index[o5*128 + 124], arr_index[o5*128 + 126], o5%2==0);
        compareSwap(arr[o5*128 + 125], arr[o5*128 + 127], arr_index[o5*128 + 125], arr_index[o5*128 + 127], o5%2==0);
    }
    for(int o6 = 0; o6 < 2; o6++){

        compareSwap(arr[o6*128], arr[o6*128 + 1], arr_index[o6*128], arr_index[o6*128 + 1], o6%2==0);
        compareSwap(arr[o6*128 + 2], arr[o6*128 + 3], arr_index[o6*128 + 2], arr_index[o6*128 + 3], o6%2==0);
        compareSwap(arr[o6*128 + 4], arr[o6*128 + 5], arr_index[o6*128 + 4], arr_index[o6*128 + 5], o6%2==0);
        compareSwap(arr[o6*128 + 6], arr[o6*128 + 7], arr_index[o6*128 + 6], arr_index[o6*128 + 7], o6%2==0);
        compareSwap(arr[o6*128 + 8], arr[o6*128 + 9], arr_index[o6*128 + 8], arr_index[o6*128 + 9], o6%2==0);
        compareSwap(arr[o6*128 + 10], arr[o6*128 + 11], arr_index[o6*128 + 10], arr_index[o6*128 + 11], o6%2==0);
        compareSwap(arr[o6*128 + 12], arr[o6*128 + 13], arr_index[o6*128 + 12], arr_index[o6*128 + 13], o6%2==0);
        compareSwap(arr[o6*128 + 14], arr[o6*128 + 15], arr_index[o6*128 + 14], arr_index[o6*128 + 15], o6%2==0);
        compareSwap(arr[o6*128 + 16], arr[o6*128 + 17], arr_index[o6*128 + 16], arr_index[o6*128 + 17], o6%2==0);
        compareSwap(arr[o6*128 + 18], arr[o6*128 + 19], arr_index[o6*128 + 18], arr_index[o6*128 + 19], o6%2==0);
        compareSwap(arr[o6*128 + 20], arr[o6*128 + 21], arr_index[o6*128 + 20], arr_index[o6*128 + 21], o6%2==0);
        compareSwap(arr[o6*128 + 22], arr[o6*128 + 23], arr_index[o6*128 + 22], arr_index[o6*128 + 23], o6%2==0);
        compareSwap(arr[o6*128 + 24], arr[o6*128 + 25], arr_index[o6*128 + 24], arr_index[o6*128 + 25], o6%2==0);
        compareSwap(arr[o6*128 + 26], arr[o6*128 + 27], arr_index[o6*128 + 26], arr_index[o6*128 + 27], o6%2==0);
        compareSwap(arr[o6*128 + 28], arr[o6*128 + 29], arr_index[o6*128 + 28], arr_index[o6*128 + 29], o6%2==0);
        compareSwap(arr[o6*128 + 30], arr[o6*128 + 31], arr_index[o6*128 + 30], arr_index[o6*128 + 31], o6%2==0);
        compareSwap(arr[o6*128 + 32], arr[o6*128 + 33], arr_index[o6*128 + 32], arr_index[o6*128 + 33], o6%2==0);
        compareSwap(arr[o6*128 + 34], arr[o6*128 + 35], arr_index[o6*128 + 34], arr_index[o6*128 + 35], o6%2==0);
        compareSwap(arr[o6*128 + 36], arr[o6*128 + 37], arr_index[o6*128 + 36], arr_index[o6*128 + 37], o6%2==0);
        compareSwap(arr[o6*128 + 38], arr[o6*128 + 39], arr_index[o6*128 + 38], arr_index[o6*128 + 39], o6%2==0);
        compareSwap(arr[o6*128 + 40], arr[o6*128 + 41], arr_index[o6*128 + 40], arr_index[o6*128 + 41], o6%2==0);
        compareSwap(arr[o6*128 + 42], arr[o6*128 + 43], arr_index[o6*128 + 42], arr_index[o6*128 + 43], o6%2==0);
        compareSwap(arr[o6*128 + 44], arr[o6*128 + 45], arr_index[o6*128 + 44], arr_index[o6*128 + 45], o6%2==0);
        compareSwap(arr[o6*128 + 46], arr[o6*128 + 47], arr_index[o6*128 + 46], arr_index[o6*128 + 47], o6%2==0);
        compareSwap(arr[o6*128 + 48], arr[o6*128 + 49], arr_index[o6*128 + 48], arr_index[o6*128 + 49], o6%2==0);
        compareSwap(arr[o6*128 + 50], arr[o6*128 + 51], arr_index[o6*128 + 50], arr_index[o6*128 + 51], o6%2==0);
        compareSwap(arr[o6*128 + 52], arr[o6*128 + 53], arr_index[o6*128 + 52], arr_index[o6*128 + 53], o6%2==0);
        compareSwap(arr[o6*128 + 54], arr[o6*128 + 55], arr_index[o6*128 + 54], arr_index[o6*128 + 55], o6%2==0);
        compareSwap(arr[o6*128 + 56], arr[o6*128 + 57], arr_index[o6*128 + 56], arr_index[o6*128 + 57], o6%2==0);
        compareSwap(arr[o6*128 + 58], arr[o6*128 + 59], arr_index[o6*128 + 58], arr_index[o6*128 + 59], o6%2==0);
        compareSwap(arr[o6*128 + 60], arr[o6*128 + 61], arr_index[o6*128 + 60], arr_index[o6*128 + 61], o6%2==0);
        compareSwap(arr[o6*128 + 62], arr[o6*128 + 63], arr_index[o6*128 + 62], arr_index[o6*128 + 63], o6%2==0);
        compareSwap(arr[o6*128 + 64], arr[o6*128 + 65], arr_index[o6*128 + 64], arr_index[o6*128 + 65], o6%2==0);
        compareSwap(arr[o6*128 + 66], arr[o6*128 + 67], arr_index[o6*128 + 66], arr_index[o6*128 + 67], o6%2==0);
        compareSwap(arr[o6*128 + 68], arr[o6*128 + 69], arr_index[o6*128 + 68], arr_index[o6*128 + 69], o6%2==0);
        compareSwap(arr[o6*128 + 70], arr[o6*128 + 71], arr_index[o6*128 + 70], arr_index[o6*128 + 71], o6%2==0);
        compareSwap(arr[o6*128 + 72], arr[o6*128 + 73], arr_index[o6*128 + 72], arr_index[o6*128 + 73], o6%2==0);
        compareSwap(arr[o6*128 + 74], arr[o6*128 + 75], arr_index[o6*128 + 74], arr_index[o6*128 + 75], o6%2==0);
        compareSwap(arr[o6*128 + 76], arr[o6*128 + 77], arr_index[o6*128 + 76], arr_index[o6*128 + 77], o6%2==0);
        compareSwap(arr[o6*128 + 78], arr[o6*128 + 79], arr_index[o6*128 + 78], arr_index[o6*128 + 79], o6%2==0);
        compareSwap(arr[o6*128 + 80], arr[o6*128 + 81], arr_index[o6*128 + 80], arr_index[o6*128 + 81], o6%2==0);
        compareSwap(arr[o6*128 + 82], arr[o6*128 + 83], arr_index[o6*128 + 82], arr_index[o6*128 + 83], o6%2==0);
        compareSwap(arr[o6*128 + 84], arr[o6*128 + 85], arr_index[o6*128 + 84], arr_index[o6*128 + 85], o6%2==0);
        compareSwap(arr[o6*128 + 86], arr[o6*128 + 87], arr_index[o6*128 + 86], arr_index[o6*128 + 87], o6%2==0);
        compareSwap(arr[o6*128 + 88], arr[o6*128 + 89], arr_index[o6*128 + 88], arr_index[o6*128 + 89], o6%2==0);
        compareSwap(arr[o6*128 + 90], arr[o6*128 + 91], arr_index[o6*128 + 90], arr_index[o6*128 + 91], o6%2==0);
        compareSwap(arr[o6*128 + 92], arr[o6*128 + 93], arr_index[o6*128 + 92], arr_index[o6*128 + 93], o6%2==0);
        compareSwap(arr[o6*128 + 94], arr[o6*128 + 95], arr_index[o6*128 + 94], arr_index[o6*128 + 95], o6%2==0);
        compareSwap(arr[o6*128 + 96], arr[o6*128 + 97], arr_index[o6*128 + 96], arr_index[o6*128 + 97], o6%2==0);
        compareSwap(arr[o6*128 + 98], arr[o6*128 + 99], arr_index[o6*128 + 98], arr_index[o6*128 + 99], o6%2==0);
        compareSwap(arr[o6*128 + 100], arr[o6*128 + 101], arr_index[o6*128 + 100], arr_index[o6*128 + 101], o6%2==0);
        compareSwap(arr[o6*128 + 102], arr[o6*128 + 103], arr_index[o6*128 + 102], arr_index[o6*128 + 103], o6%2==0);
        compareSwap(arr[o6*128 + 104], arr[o6*128 + 105], arr_index[o6*128 + 104], arr_index[o6*128 + 105], o6%2==0);
        compareSwap(arr[o6*128 + 106], arr[o6*128 + 107], arr_index[o6*128 + 106], arr_index[o6*128 + 107], o6%2==0);
        compareSwap(arr[o6*128 + 108], arr[o6*128 + 109], arr_index[o6*128 + 108], arr_index[o6*128 + 109], o6%2==0);
        compareSwap(arr[o6*128 + 110], arr[o6*128 + 111], arr_index[o6*128 + 110], arr_index[o6*128 + 111], o6%2==0);
        compareSwap(arr[o6*128 + 112], arr[o6*128 + 113], arr_index[o6*128 + 112], arr_index[o6*128 + 113], o6%2==0);
        compareSwap(arr[o6*128 + 114], arr[o6*128 + 115], arr_index[o6*128 + 114], arr_index[o6*128 + 115], o6%2==0);
        compareSwap(arr[o6*128 + 116], arr[o6*128 + 117], arr_index[o6*128 + 116], arr_index[o6*128 + 117], o6%2==0);
        compareSwap(arr[o6*128 + 118], arr[o6*128 + 119], arr_index[o6*128 + 118], arr_index[o6*128 + 119], o6%2==0);
        compareSwap(arr[o6*128 + 120], arr[o6*128 + 121], arr_index[o6*128 + 120], arr_index[o6*128 + 121], o6%2==0);
        compareSwap(arr[o6*128 + 122], arr[o6*128 + 123], arr_index[o6*128 + 122], arr_index[o6*128 + 123], o6%2==0);
        compareSwap(arr[o6*128 + 124], arr[o6*128 + 125], arr_index[o6*128 + 124], arr_index[o6*128 + 125], o6%2==0);
        compareSwap(arr[o6*128 + 126], arr[o6*128 + 127], arr_index[o6*128 + 126], arr_index[o6*128 + 127], o6%2==0);
    }
    //stage 8
    for(int p = 0; p < 1; p++){
        compareSwap(arr[p*256 + 0], arr[p*256 + 128], arr_index[p*256 + 0], arr_index[p*256 + 128], p%2==0);
        compareSwap(arr[p*256 + 1], arr[p*256 + 129], arr_index[p*256 + 1], arr_index[p*256 + 129], p%2==0);
        compareSwap(arr[p*256 + 2], arr[p*256 + 130], arr_index[p*256 + 2], arr_index[p*256 + 130], p%2==0);
        compareSwap(arr[p*256 + 3], arr[p*256 + 131], arr_index[p*256 + 3], arr_index[p*256 + 131], p%2==0);
        compareSwap(arr[p*256 + 4], arr[p*256 + 132], arr_index[p*256 + 4], arr_index[p*256 + 132], p%2==0);
        compareSwap(arr[p*256 + 5], arr[p*256 + 133], arr_index[p*256 + 5], arr_index[p*256 + 133], p%2==0);
        compareSwap(arr[p*256 + 6], arr[p*256 + 134], arr_index[p*256 + 6], arr_index[p*256 + 134], p%2==0);
        compareSwap(arr[p*256 + 7], arr[p*256 + 135], arr_index[p*256 + 7], arr_index[p*256 + 135], p%2==0);
        compareSwap(arr[p*256 + 8], arr[p*256 + 136], arr_index[p*256 + 8], arr_index[p*256 + 136], p%2==0);
        compareSwap(arr[p*256 + 9], arr[p*256 + 137], arr_index[p*256 + 9], arr_index[p*256 + 137], p%2==0);
        compareSwap(arr[p*256 + 10], arr[p*256 + 138], arr_index[p*256 + 10], arr_index[p*256 + 138], p%2==0);
        compareSwap(arr[p*256 + 11], arr[p*256 + 139], arr_index[p*256 + 11], arr_index[p*256 + 139], p%2==0);
        compareSwap(arr[p*256 + 12], arr[p*256 + 140], arr_index[p*256 + 12], arr_index[p*256 + 140], p%2==0);
        compareSwap(arr[p*256 + 13], arr[p*256 + 141], arr_index[p*256 + 13], arr_index[p*256 + 141], p%2==0);
        compareSwap(arr[p*256 + 14], arr[p*256 + 142], arr_index[p*256 + 14], arr_index[p*256 + 142], p%2==0);
        compareSwap(arr[p*256 + 15], arr[p*256 + 143], arr_index[p*256 + 15], arr_index[p*256 + 143], p%2==0);
        compareSwap(arr[p*256 + 16], arr[p*256 + 144], arr_index[p*256 + 16], arr_index[p*256 + 144], p%2==0);
        compareSwap(arr[p*256 + 17], arr[p*256 + 145], arr_index[p*256 + 17], arr_index[p*256 + 145], p%2==0);
        compareSwap(arr[p*256 + 18], arr[p*256 + 146], arr_index[p*256 + 18], arr_index[p*256 + 146], p%2==0);
        compareSwap(arr[p*256 + 19], arr[p*256 + 147], arr_index[p*256 + 19], arr_index[p*256 + 147], p%2==0);
        compareSwap(arr[p*256 + 20], arr[p*256 + 148], arr_index[p*256 + 20], arr_index[p*256 + 148], p%2==0);
        compareSwap(arr[p*256 + 21], arr[p*256 + 149], arr_index[p*256 + 21], arr_index[p*256 + 149], p%2==0);
        compareSwap(arr[p*256 + 22], arr[p*256 + 150], arr_index[p*256 + 22], arr_index[p*256 + 150], p%2==0);
        compareSwap(arr[p*256 + 23], arr[p*256 + 151], arr_index[p*256 + 23], arr_index[p*256 + 151], p%2==0);
        compareSwap(arr[p*256 + 24], arr[p*256 + 152], arr_index[p*256 + 24], arr_index[p*256 + 152], p%2==0);
        compareSwap(arr[p*256 + 25], arr[p*256 + 153], arr_index[p*256 + 25], arr_index[p*256 + 153], p%2==0);
        compareSwap(arr[p*256 + 26], arr[p*256 + 154], arr_index[p*256 + 26], arr_index[p*256 + 154], p%2==0);
        compareSwap(arr[p*256 + 27], arr[p*256 + 155], arr_index[p*256 + 27], arr_index[p*256 + 155], p%2==0);
        compareSwap(arr[p*256 + 28], arr[p*256 + 156], arr_index[p*256 + 28], arr_index[p*256 + 156], p%2==0);
        compareSwap(arr[p*256 + 29], arr[p*256 + 157], arr_index[p*256 + 29], arr_index[p*256 + 157], p%2==0);
        compareSwap(arr[p*256 + 30], arr[p*256 + 158], arr_index[p*256 + 30], arr_index[p*256 + 158], p%2==0);
        compareSwap(arr[p*256 + 31], arr[p*256 + 159], arr_index[p*256 + 31], arr_index[p*256 + 159], p%2==0);
        compareSwap(arr[p*256 + 32], arr[p*256 + 160], arr_index[p*256 + 32], arr_index[p*256 + 160], p%2==0);
        compareSwap(arr[p*256 + 33], arr[p*256 + 161], arr_index[p*256 + 33], arr_index[p*256 + 161], p%2==0);
        compareSwap(arr[p*256 + 34], arr[p*256 + 162], arr_index[p*256 + 34], arr_index[p*256 + 162], p%2==0);
        compareSwap(arr[p*256 + 35], arr[p*256 + 163], arr_index[p*256 + 35], arr_index[p*256 + 163], p%2==0);
        compareSwap(arr[p*256 + 36], arr[p*256 + 164], arr_index[p*256 + 36], arr_index[p*256 + 164], p%2==0);
        compareSwap(arr[p*256 + 37], arr[p*256 + 165], arr_index[p*256 + 37], arr_index[p*256 + 165], p%2==0);
        compareSwap(arr[p*256 + 38], arr[p*256 + 166], arr_index[p*256 + 38], arr_index[p*256 + 166], p%2==0);
        compareSwap(arr[p*256 + 39], arr[p*256 + 167], arr_index[p*256 + 39], arr_index[p*256 + 167], p%2==0);
        compareSwap(arr[p*256 + 40], arr[p*256 + 168], arr_index[p*256 + 40], arr_index[p*256 + 168], p%2==0);
        compareSwap(arr[p*256 + 41], arr[p*256 + 169], arr_index[p*256 + 41], arr_index[p*256 + 169], p%2==0);
        compareSwap(arr[p*256 + 42], arr[p*256 + 170], arr_index[p*256 + 42], arr_index[p*256 + 170], p%2==0);
        compareSwap(arr[p*256 + 43], arr[p*256 + 171], arr_index[p*256 + 43], arr_index[p*256 + 171], p%2==0);
        compareSwap(arr[p*256 + 44], arr[p*256 + 172], arr_index[p*256 + 44], arr_index[p*256 + 172], p%2==0);
        compareSwap(arr[p*256 + 45], arr[p*256 + 173], arr_index[p*256 + 45], arr_index[p*256 + 173], p%2==0);
        compareSwap(arr[p*256 + 46], arr[p*256 + 174], arr_index[p*256 + 46], arr_index[p*256 + 174], p%2==0);
        compareSwap(arr[p*256 + 47], arr[p*256 + 175], arr_index[p*256 + 47], arr_index[p*256 + 175], p%2==0);
        compareSwap(arr[p*256 + 48], arr[p*256 + 176], arr_index[p*256 + 48], arr_index[p*256 + 176], p%2==0);
        compareSwap(arr[p*256 + 49], arr[p*256 + 177], arr_index[p*256 + 49], arr_index[p*256 + 177], p%2==0);
        compareSwap(arr[p*256 + 50], arr[p*256 + 178], arr_index[p*256 + 50], arr_index[p*256 + 178], p%2==0);
        compareSwap(arr[p*256 + 51], arr[p*256 + 179], arr_index[p*256 + 51], arr_index[p*256 + 179], p%2==0);
        compareSwap(arr[p*256 + 52], arr[p*256 + 180], arr_index[p*256 + 52], arr_index[p*256 + 180], p%2==0);
        compareSwap(arr[p*256 + 53], arr[p*256 + 181], arr_index[p*256 + 53], arr_index[p*256 + 181], p%2==0);
        compareSwap(arr[p*256 + 54], arr[p*256 + 182], arr_index[p*256 + 54], arr_index[p*256 + 182], p%2==0);
        compareSwap(arr[p*256 + 55], arr[p*256 + 183], arr_index[p*256 + 55], arr_index[p*256 + 183], p%2==0);
        compareSwap(arr[p*256 + 56], arr[p*256 + 184], arr_index[p*256 + 56], arr_index[p*256 + 184], p%2==0);
        compareSwap(arr[p*256 + 57], arr[p*256 + 185], arr_index[p*256 + 57], arr_index[p*256 + 185], p%2==0);
        compareSwap(arr[p*256 + 58], arr[p*256 + 186], arr_index[p*256 + 58], arr_index[p*256 + 186], p%2==0);
        compareSwap(arr[p*256 + 59], arr[p*256 + 187], arr_index[p*256 + 59], arr_index[p*256 + 187], p%2==0);
        compareSwap(arr[p*256 + 60], arr[p*256 + 188], arr_index[p*256 + 60], arr_index[p*256 + 188], p%2==0);
        compareSwap(arr[p*256 + 61], arr[p*256 + 189], arr_index[p*256 + 61], arr_index[p*256 + 189], p%2==0);
        compareSwap(arr[p*256 + 62], arr[p*256 + 190], arr_index[p*256 + 62], arr_index[p*256 + 190], p%2==0);
        compareSwap(arr[p*256 + 63], arr[p*256 + 191], arr_index[p*256 + 63], arr_index[p*256 + 191], p%2==0);
        compareSwap(arr[p*256 + 64], arr[p*256 + 192], arr_index[p*256 + 64], arr_index[p*256 + 192], p%2==0);
        compareSwap(arr[p*256 + 65], arr[p*256 + 193], arr_index[p*256 + 65], arr_index[p*256 + 193], p%2==0);
        compareSwap(arr[p*256 + 66], arr[p*256 + 194], arr_index[p*256 + 66], arr_index[p*256 + 194], p%2==0);
        compareSwap(arr[p*256 + 67], arr[p*256 + 195], arr_index[p*256 + 67], arr_index[p*256 + 195], p%2==0);
        compareSwap(arr[p*256 + 68], arr[p*256 + 196], arr_index[p*256 + 68], arr_index[p*256 + 196], p%2==0);
        compareSwap(arr[p*256 + 69], arr[p*256 + 197], arr_index[p*256 + 69], arr_index[p*256 + 197], p%2==0);
        compareSwap(arr[p*256 + 70], arr[p*256 + 198], arr_index[p*256 + 70], arr_index[p*256 + 198], p%2==0);
        compareSwap(arr[p*256 + 71], arr[p*256 + 199], arr_index[p*256 + 71], arr_index[p*256 + 199], p%2==0);
        compareSwap(arr[p*256 + 72], arr[p*256 + 200], arr_index[p*256 + 72], arr_index[p*256 + 200], p%2==0);
        compareSwap(arr[p*256 + 73], arr[p*256 + 201], arr_index[p*256 + 73], arr_index[p*256 + 201], p%2==0);
        compareSwap(arr[p*256 + 74], arr[p*256 + 202], arr_index[p*256 + 74], arr_index[p*256 + 202], p%2==0);
        compareSwap(arr[p*256 + 75], arr[p*256 + 203], arr_index[p*256 + 75], arr_index[p*256 + 203], p%2==0);
        compareSwap(arr[p*256 + 76], arr[p*256 + 204], arr_index[p*256 + 76], arr_index[p*256 + 204], p%2==0);
        compareSwap(arr[p*256 + 77], arr[p*256 + 205], arr_index[p*256 + 77], arr_index[p*256 + 205], p%2==0);
        compareSwap(arr[p*256 + 78], arr[p*256 + 206], arr_index[p*256 + 78], arr_index[p*256 + 206], p%2==0);
        compareSwap(arr[p*256 + 79], arr[p*256 + 207], arr_index[p*256 + 79], arr_index[p*256 + 207], p%2==0);
        compareSwap(arr[p*256 + 80], arr[p*256 + 208], arr_index[p*256 + 80], arr_index[p*256 + 208], p%2==0);
        compareSwap(arr[p*256 + 81], arr[p*256 + 209], arr_index[p*256 + 81], arr_index[p*256 + 209], p%2==0);
        compareSwap(arr[p*256 + 82], arr[p*256 + 210], arr_index[p*256 + 82], arr_index[p*256 + 210], p%2==0);
        compareSwap(arr[p*256 + 83], arr[p*256 + 211], arr_index[p*256 + 83], arr_index[p*256 + 211], p%2==0);
        compareSwap(arr[p*256 + 84], arr[p*256 + 212], arr_index[p*256 + 84], arr_index[p*256 + 212], p%2==0);
        compareSwap(arr[p*256 + 85], arr[p*256 + 213], arr_index[p*256 + 85], arr_index[p*256 + 213], p%2==0);
        compareSwap(arr[p*256 + 86], arr[p*256 + 214], arr_index[p*256 + 86], arr_index[p*256 + 214], p%2==0);
        compareSwap(arr[p*256 + 87], arr[p*256 + 215], arr_index[p*256 + 87], arr_index[p*256 + 215], p%2==0);
        compareSwap(arr[p*256 + 88], arr[p*256 + 216], arr_index[p*256 + 88], arr_index[p*256 + 216], p%2==0);
        compareSwap(arr[p*256 + 89], arr[p*256 + 217], arr_index[p*256 + 89], arr_index[p*256 + 217], p%2==0);
        compareSwap(arr[p*256 + 90], arr[p*256 + 218], arr_index[p*256 + 90], arr_index[p*256 + 218], p%2==0);
        compareSwap(arr[p*256 + 91], arr[p*256 + 219], arr_index[p*256 + 91], arr_index[p*256 + 219], p%2==0);
        compareSwap(arr[p*256 + 92], arr[p*256 + 220], arr_index[p*256 + 92], arr_index[p*256 + 220], p%2==0);
        compareSwap(arr[p*256 + 93], arr[p*256 + 221], arr_index[p*256 + 93], arr_index[p*256 + 221], p%2==0);
        compareSwap(arr[p*256 + 94], arr[p*256 + 222], arr_index[p*256 + 94], arr_index[p*256 + 222], p%2==0);
        compareSwap(arr[p*256 + 95], arr[p*256 + 223], arr_index[p*256 + 95], arr_index[p*256 + 223], p%2==0);
        compareSwap(arr[p*256 + 96], arr[p*256 + 224], arr_index[p*256 + 96], arr_index[p*256 + 224], p%2==0);
        compareSwap(arr[p*256 + 97], arr[p*256 + 225], arr_index[p*256 + 97], arr_index[p*256 + 225], p%2==0);
        compareSwap(arr[p*256 + 98], arr[p*256 + 226], arr_index[p*256 + 98], arr_index[p*256 + 226], p%2==0);
        compareSwap(arr[p*256 + 99], arr[p*256 + 227], arr_index[p*256 + 99], arr_index[p*256 + 227], p%2==0);
        compareSwap(arr[p*256 + 100], arr[p*256 + 228], arr_index[p*256 + 100], arr_index[p*256 + 228], p%2==0);
        compareSwap(arr[p*256 + 101], arr[p*256 + 229], arr_index[p*256 + 101], arr_index[p*256 + 229], p%2==0);
        compareSwap(arr[p*256 + 102], arr[p*256 + 230], arr_index[p*256 + 102], arr_index[p*256 + 230], p%2==0);
        compareSwap(arr[p*256 + 103], arr[p*256 + 231], arr_index[p*256 + 103], arr_index[p*256 + 231], p%2==0);
        compareSwap(arr[p*256 + 104], arr[p*256 + 232], arr_index[p*256 + 104], arr_index[p*256 + 232], p%2==0);
        compareSwap(arr[p*256 + 105], arr[p*256 + 233], arr_index[p*256 + 105], arr_index[p*256 + 233], p%2==0);
        compareSwap(arr[p*256 + 106], arr[p*256 + 234], arr_index[p*256 + 106], arr_index[p*256 + 234], p%2==0);
        compareSwap(arr[p*256 + 107], arr[p*256 + 235], arr_index[p*256 + 107], arr_index[p*256 + 235], p%2==0);
        compareSwap(arr[p*256 + 108], arr[p*256 + 236], arr_index[p*256 + 108], arr_index[p*256 + 236], p%2==0);
        compareSwap(arr[p*256 + 109], arr[p*256 + 237], arr_index[p*256 + 109], arr_index[p*256 + 237], p%2==0);
        compareSwap(arr[p*256 + 110], arr[p*256 + 238], arr_index[p*256 + 110], arr_index[p*256 + 238], p%2==0);
        compareSwap(arr[p*256 + 111], arr[p*256 + 239], arr_index[p*256 + 111], arr_index[p*256 + 239], p%2==0);
        compareSwap(arr[p*256 + 112], arr[p*256 + 240], arr_index[p*256 + 112], arr_index[p*256 + 240], p%2==0);
        compareSwap(arr[p*256 + 113], arr[p*256 + 241], arr_index[p*256 + 113], arr_index[p*256 + 241], p%2==0);
        compareSwap(arr[p*256 + 114], arr[p*256 + 242], arr_index[p*256 + 114], arr_index[p*256 + 242], p%2==0);
        compareSwap(arr[p*256 + 115], arr[p*256 + 243], arr_index[p*256 + 115], arr_index[p*256 + 243], p%2==0);
        compareSwap(arr[p*256 + 116], arr[p*256 + 244], arr_index[p*256 + 116], arr_index[p*256 + 244], p%2==0);
        compareSwap(arr[p*256 + 117], arr[p*256 + 245], arr_index[p*256 + 117], arr_index[p*256 + 245], p%2==0);
        compareSwap(arr[p*256 + 118], arr[p*256 + 246], arr_index[p*256 + 118], arr_index[p*256 + 246], p%2==0);
        compareSwap(arr[p*256 + 119], arr[p*256 + 247], arr_index[p*256 + 119], arr_index[p*256 + 247], p%2==0);
        compareSwap(arr[p*256 + 120], arr[p*256 + 248], arr_index[p*256 + 120], arr_index[p*256 + 248], p%2==0);
        compareSwap(arr[p*256 + 121], arr[p*256 + 249], arr_index[p*256 + 121], arr_index[p*256 + 249], p%2==0);
        compareSwap(arr[p*256 + 122], arr[p*256 + 250], arr_index[p*256 + 122], arr_index[p*256 + 250], p%2==0);
        compareSwap(arr[p*256 + 123], arr[p*256 + 251], arr_index[p*256 + 123], arr_index[p*256 + 251], p%2==0);
        compareSwap(arr[p*256 + 124], arr[p*256 + 252], arr_index[p*256 + 124], arr_index[p*256 + 252], p%2==0);
        compareSwap(arr[p*256 + 125], arr[p*256 + 253], arr_index[p*256 + 125], arr_index[p*256 + 253], p%2==0);
        compareSwap(arr[p*256 + 126], arr[p*256 + 254], arr_index[p*256 + 126], arr_index[p*256 + 254], p%2==0);
        compareSwap(arr[p*256 + 127], arr[p*256 + 255], arr_index[p*256 + 127], arr_index[p*256 + 255], p%2==0);
    }

    for(int p1 = 0; p1 < 1; p1++)
    {
        compareSwap(arr[p1*256 + 0], arr[p1*256 + 64], arr_index[p1*256 + 0], arr_index[p1*256 + 64], p1%2==0);
        compareSwap(arr[p1*256 + 1], arr[p1*256 + 65], arr_index[p1*256 + 1], arr_index[p1*256 + 65], p1%2==0);
        compareSwap(arr[p1*256 + 2], arr[p1*256 + 66], arr_index[p1*256 + 2], arr_index[p1*256 + 66], p1%2==0);
        compareSwap(arr[p1*256 + 3], arr[p1*256 + 67], arr_index[p1*256 + 3], arr_index[p1*256 + 67], p1%2==0);
        compareSwap(arr[p1*256 + 4], arr[p1*256 + 68], arr_index[p1*256 + 4], arr_index[p1*256 + 68], p1%2==0);
        compareSwap(arr[p1*256 + 5], arr[p1*256 + 69], arr_index[p1*256 + 5], arr_index[p1*256 + 69], p1%2==0);
        compareSwap(arr[p1*256 + 6], arr[p1*256 + 70], arr_index[p1*256 + 6], arr_index[p1*256 + 70], p1%2==0);
        compareSwap(arr[p1*256 + 7], arr[p1*256 + 71], arr_index[p1*256 + 7], arr_index[p1*256 + 71], p1%2==0);
        compareSwap(arr[p1*256 + 8], arr[p1*256 + 72], arr_index[p1*256 + 8], arr_index[p1*256 + 72], p1%2==0);
        compareSwap(arr[p1*256 + 9], arr[p1*256 + 73], arr_index[p1*256 + 9], arr_index[p1*256 + 73], p1%2==0);
        compareSwap(arr[p1*256 + 10], arr[p1*256 + 74], arr_index[p1*256 + 10], arr_index[p1*256 + 74], p1%2==0);
        compareSwap(arr[p1*256 + 11], arr[p1*256 + 75], arr_index[p1*256 + 11], arr_index[p1*256 + 75], p1%2==0);
        compareSwap(arr[p1*256 + 12], arr[p1*256 + 76], arr_index[p1*256 + 12], arr_index[p1*256 + 76], p1%2==0);
        compareSwap(arr[p1*256 + 13], arr[p1*256 + 77], arr_index[p1*256 + 13], arr_index[p1*256 + 77], p1%2==0);
        compareSwap(arr[p1*256 + 14], arr[p1*256 + 78], arr_index[p1*256 + 14], arr_index[p1*256 + 78], p1%2==0);
        compareSwap(arr[p1*256 + 15], arr[p1*256 + 79], arr_index[p1*256 + 15], arr_index[p1*256 + 79], p1%2==0);
        compareSwap(arr[p1*256 + 16], arr[p1*256 + 80], arr_index[p1*256 + 16], arr_index[p1*256 + 80], p1%2==0);
        compareSwap(arr[p1*256 + 17], arr[p1*256 + 81], arr_index[p1*256 + 17], arr_index[p1*256 + 81], p1%2==0);
        compareSwap(arr[p1*256 + 18], arr[p1*256 + 82], arr_index[p1*256 + 18], arr_index[p1*256 + 82], p1%2==0);
        compareSwap(arr[p1*256 + 19], arr[p1*256 + 83], arr_index[p1*256 + 19], arr_index[p1*256 + 83], p1%2==0);
        compareSwap(arr[p1*256 + 20], arr[p1*256 + 84], arr_index[p1*256 + 20], arr_index[p1*256 + 84], p1%2==0);
        compareSwap(arr[p1*256 + 21], arr[p1*256 + 85], arr_index[p1*256 + 21], arr_index[p1*256 + 85], p1%2==0);
        compareSwap(arr[p1*256 + 22], arr[p1*256 + 86], arr_index[p1*256 + 22], arr_index[p1*256 + 86], p1%2==0);
        compareSwap(arr[p1*256 + 23], arr[p1*256 + 87], arr_index[p1*256 + 23], arr_index[p1*256 + 87], p1%2==0);
        compareSwap(arr[p1*256 + 24], arr[p1*256 + 88], arr_index[p1*256 + 24], arr_index[p1*256 + 88], p1%2==0);
        compareSwap(arr[p1*256 + 25], arr[p1*256 + 89], arr_index[p1*256 + 25], arr_index[p1*256 + 89], p1%2==0);
        compareSwap(arr[p1*256 + 26], arr[p1*256 + 90], arr_index[p1*256 + 26], arr_index[p1*256 + 90], p1%2==0);
        compareSwap(arr[p1*256 + 27], arr[p1*256 + 91], arr_index[p1*256 + 27], arr_index[p1*256 + 91], p1%2==0);
        compareSwap(arr[p1*256 + 28], arr[p1*256 + 92], arr_index[p1*256 + 28], arr_index[p1*256 + 92], p1%2==0);
        compareSwap(arr[p1*256 + 29], arr[p1*256 + 93], arr_index[p1*256 + 29], arr_index[p1*256 + 93], p1%2==0);
        compareSwap(arr[p1*256 + 30], arr[p1*256 + 94], arr_index[p1*256 + 30], arr_index[p1*256 + 94], p1%2==0);
        compareSwap(arr[p1*256 + 31], arr[p1*256 + 95], arr_index[p1*256 + 31], arr_index[p1*256 + 95], p1%2==0);
        compareSwap(arr[p1*256 + 32], arr[p1*256 + 96], arr_index[p1*256 + 32], arr_index[p1*256 + 96], p1%2==0);
        compareSwap(arr[p1*256 + 33], arr[p1*256 + 97], arr_index[p1*256 + 33], arr_index[p1*256 + 97], p1%2==0);
        compareSwap(arr[p1*256 + 34], arr[p1*256 + 98], arr_index[p1*256 + 34], arr_index[p1*256 + 98], p1%2==0);
        compareSwap(arr[p1*256 + 35], arr[p1*256 + 99], arr_index[p1*256 + 35], arr_index[p1*256 + 99], p1%2==0);
        compareSwap(arr[p1*256 + 36], arr[p1*256 + 100], arr_index[p1*256 + 36], arr_index[p1*256 + 100], p1%2==0);
        compareSwap(arr[p1*256 + 37], arr[p1*256 + 101], arr_index[p1*256 + 37], arr_index[p1*256 + 101], p1%2==0);
        compareSwap(arr[p1*256 + 38], arr[p1*256 + 102], arr_index[p1*256 + 38], arr_index[p1*256 + 102], p1%2==0);
        compareSwap(arr[p1*256 + 39], arr[p1*256 + 103], arr_index[p1*256 + 39], arr_index[p1*256 + 103], p1%2==0);
        compareSwap(arr[p1*256 + 40], arr[p1*256 + 104], arr_index[p1*256 + 40], arr_index[p1*256 + 104], p1%2==0);
        compareSwap(arr[p1*256 + 41], arr[p1*256 + 105], arr_index[p1*256 + 41], arr_index[p1*256 + 105], p1%2==0);
        compareSwap(arr[p1*256 + 42], arr[p1*256 + 106], arr_index[p1*256 + 42], arr_index[p1*256 + 106], p1%2==0);
        compareSwap(arr[p1*256 + 43], arr[p1*256 + 107], arr_index[p1*256 + 43], arr_index[p1*256 + 107], p1%2==0);
        compareSwap(arr[p1*256 + 44], arr[p1*256 + 108], arr_index[p1*256 + 44], arr_index[p1*256 + 108], p1%2==0);
        compareSwap(arr[p1*256 + 45], arr[p1*256 + 109], arr_index[p1*256 + 45], arr_index[p1*256 + 109], p1%2==0);
        compareSwap(arr[p1*256 + 46], arr[p1*256 + 110], arr_index[p1*256 + 46], arr_index[p1*256 + 110], p1%2==0);
        compareSwap(arr[p1*256 + 47], arr[p1*256 + 111], arr_index[p1*256 + 47], arr_index[p1*256 + 111], p1%2==0);
        compareSwap(arr[p1*256 + 48], arr[p1*256 + 112], arr_index[p1*256 + 48], arr_index[p1*256 + 112], p1%2==0);
        compareSwap(arr[p1*256 + 49], arr[p1*256 + 113], arr_index[p1*256 + 49], arr_index[p1*256 + 113], p1%2==0);
        compareSwap(arr[p1*256 + 50], arr[p1*256 + 114], arr_index[p1*256 + 50], arr_index[p1*256 + 114], p1%2==0);
        compareSwap(arr[p1*256 + 51], arr[p1*256 + 115], arr_index[p1*256 + 51], arr_index[p1*256 + 115], p1%2==0);
        compareSwap(arr[p1*256 + 52], arr[p1*256 + 116], arr_index[p1*256 + 52], arr_index[p1*256 + 116], p1%2==0);
        compareSwap(arr[p1*256 + 53], arr[p1*256 + 117], arr_index[p1*256 + 53], arr_index[p1*256 + 117], p1%2==0);
        compareSwap(arr[p1*256 + 54], arr[p1*256 + 118], arr_index[p1*256 + 54], arr_index[p1*256 + 118], p1%2==0);
        compareSwap(arr[p1*256 + 55], arr[p1*256 + 119], arr_index[p1*256 + 55], arr_index[p1*256 + 119], p1%2==0);
        compareSwap(arr[p1*256 + 56], arr[p1*256 + 120], arr_index[p1*256 + 56], arr_index[p1*256 + 120], p1%2==0);
        compareSwap(arr[p1*256 + 57], arr[p1*256 + 121], arr_index[p1*256 + 57], arr_index[p1*256 + 121], p1%2==0);
        compareSwap(arr[p1*256 + 58], arr[p1*256 + 122], arr_index[p1*256 + 58], arr_index[p1*256 + 122], p1%2==0);
        compareSwap(arr[p1*256 + 59], arr[p1*256 + 123], arr_index[p1*256 + 59], arr_index[p1*256 + 123], p1%2==0);
        compareSwap(arr[p1*256 + 60], arr[p1*256 + 124], arr_index[p1*256 + 60], arr_index[p1*256 + 124], p1%2==0);
        compareSwap(arr[p1*256 + 61], arr[p1*256 + 125], arr_index[p1*256 + 61], arr_index[p1*256 + 125], p1%2==0);
        compareSwap(arr[p1*256 + 62], arr[p1*256 + 126], arr_index[p1*256 + 62], arr_index[p1*256 + 126], p1%2==0);
        compareSwap(arr[p1*256 + 63], arr[p1*256 + 127], arr_index[p1*256 + 63], arr_index[p1*256 + 127], p1%2==0);

        compareSwap(arr[p1*256 + 128], arr[p1*256 + 192], arr_index[p1*256 + 128], arr_index[p1*256 + 192], p1%2==0);
        compareSwap(arr[p1*256 + 129], arr[p1*256 + 193], arr_index[p1*256 + 129], arr_index[p1*256 + 193], p1%2==0);
        compareSwap(arr[p1*256 + 130], arr[p1*256 + 194], arr_index[p1*256 + 130], arr_index[p1*256 + 194], p1%2==0);
        compareSwap(arr[p1*256 + 131], arr[p1*256 + 195], arr_index[p1*256 + 131], arr_index[p1*256 + 195], p1%2==0);
        compareSwap(arr[p1*256 + 132], arr[p1*256 + 196], arr_index[p1*256 + 132], arr_index[p1*256 + 196], p1%2==0);
        compareSwap(arr[p1*256 + 133], arr[p1*256 + 197], arr_index[p1*256 + 133], arr_index[p1*256 + 197], p1%2==0);
        compareSwap(arr[p1*256 + 134], arr[p1*256 + 198], arr_index[p1*256 + 134], arr_index[p1*256 + 198], p1%2==0);
        compareSwap(arr[p1*256 + 135], arr[p1*256 + 199], arr_index[p1*256 + 135], arr_index[p1*256 + 199], p1%2==0);
        compareSwap(arr[p1*256 + 136], arr[p1*256 + 200], arr_index[p1*256 + 136], arr_index[p1*256 + 200], p1%2==0);
        compareSwap(arr[p1*256 + 137], arr[p1*256 + 201], arr_index[p1*256 + 137], arr_index[p1*256 + 201], p1%2==0);
        compareSwap(arr[p1*256 + 138], arr[p1*256 + 202], arr_index[p1*256 + 138], arr_index[p1*256 + 202], p1%2==0);
        compareSwap(arr[p1*256 + 139], arr[p1*256 + 203], arr_index[p1*256 + 139], arr_index[p1*256 + 203], p1%2==0);
        compareSwap(arr[p1*256 + 140], arr[p1*256 + 204], arr_index[p1*256 + 140], arr_index[p1*256 + 204], p1%2==0);
        compareSwap(arr[p1*256 + 141], arr[p1*256 + 205], arr_index[p1*256 + 141], arr_index[p1*256 + 205], p1%2==0);
        compareSwap(arr[p1*256 + 142], arr[p1*256 + 206], arr_index[p1*256 + 142], arr_index[p1*256 + 206], p1%2==0);
        compareSwap(arr[p1*256 + 143], arr[p1*256 + 207], arr_index[p1*256 + 143], arr_index[p1*256 + 207], p1%2==0);
        compareSwap(arr[p1*256 + 144], arr[p1*256 + 208], arr_index[p1*256 + 144], arr_index[p1*256 + 208], p1%2==0);
        compareSwap(arr[p1*256 + 145], arr[p1*256 + 209], arr_index[p1*256 + 145], arr_index[p1*256 + 209], p1%2==0);
        compareSwap(arr[p1*256 + 146], arr[p1*256 + 210], arr_index[p1*256 + 146], arr_index[p1*256 + 210], p1%2==0);
        compareSwap(arr[p1*256 + 147], arr[p1*256 + 211], arr_index[p1*256 + 147], arr_index[p1*256 + 211], p1%2==0);
        compareSwap(arr[p1*256 + 148], arr[p1*256 + 212], arr_index[p1*256 + 148], arr_index[p1*256 + 212], p1%2==0);
        compareSwap(arr[p1*256 + 149], arr[p1*256 + 213], arr_index[p1*256 + 149], arr_index[p1*256 + 213], p1%2==0);
        compareSwap(arr[p1*256 + 150], arr[p1*256 + 214], arr_index[p1*256 + 150], arr_index[p1*256 + 214], p1%2==0);
        compareSwap(arr[p1*256 + 151], arr[p1*256 + 215], arr_index[p1*256 + 151], arr_index[p1*256 + 215], p1%2==0);
        compareSwap(arr[p1*256 + 152], arr[p1*256 + 216], arr_index[p1*256 + 152], arr_index[p1*256 + 216], p1%2==0);
        compareSwap(arr[p1*256 + 153], arr[p1*256 + 217], arr_index[p1*256 + 153], arr_index[p1*256 + 217], p1%2==0);
        compareSwap(arr[p1*256 + 154], arr[p1*256 + 218], arr_index[p1*256 + 154], arr_index[p1*256 + 218], p1%2==0);
        compareSwap(arr[p1*256 + 155], arr[p1*256 + 219], arr_index[p1*256 + 155], arr_index[p1*256 + 219], p1%2==0);
        compareSwap(arr[p1*256 + 156], arr[p1*256 + 220], arr_index[p1*256 + 156], arr_index[p1*256 + 220], p1%2==0);
        compareSwap(arr[p1*256 + 157], arr[p1*256 + 221], arr_index[p1*256 + 157], arr_index[p1*256 + 221], p1%2==0);
        compareSwap(arr[p1*256 + 158], arr[p1*256 + 222], arr_index[p1*256 + 158], arr_index[p1*256 + 222], p1%2==0);
        compareSwap(arr[p1*256 + 159], arr[p1*256 + 223], arr_index[p1*256 + 159], arr_index[p1*256 + 223], p1%2==0);
        compareSwap(arr[p1*256 + 160], arr[p1*256 + 224], arr_index[p1*256 + 160], arr_index[p1*256 + 224], p1%2==0);
        compareSwap(arr[p1*256 + 161], arr[p1*256 + 225], arr_index[p1*256 + 161], arr_index[p1*256 + 225], p1%2==0);
        compareSwap(arr[p1*256 + 162], arr[p1*256 + 226], arr_index[p1*256 + 162], arr_index[p1*256 + 226], p1%2==0);
        compareSwap(arr[p1*256 + 163], arr[p1*256 + 227], arr_index[p1*256 + 163], arr_index[p1*256 + 227], p1%2==0);
        compareSwap(arr[p1*256 + 164], arr[p1*256 + 228], arr_index[p1*256 + 164], arr_index[p1*256 + 228], p1%2==0);
        compareSwap(arr[p1*256 + 165], arr[p1*256 + 229], arr_index[p1*256 + 165], arr_index[p1*256 + 229], p1%2==0);
        compareSwap(arr[p1*256 + 166], arr[p1*256 + 230], arr_index[p1*256 + 166], arr_index[p1*256 + 230], p1%2==0);
        compareSwap(arr[p1*256 + 167], arr[p1*256 + 231], arr_index[p1*256 + 167], arr_index[p1*256 + 231], p1%2==0);
        compareSwap(arr[p1*256 + 168], arr[p1*256 + 232], arr_index[p1*256 + 168], arr_index[p1*256 + 232], p1%2==0);
        compareSwap(arr[p1*256 + 169], arr[p1*256 + 233], arr_index[p1*256 + 169], arr_index[p1*256 + 233], p1%2==0);
        compareSwap(arr[p1*256 + 170], arr[p1*256 + 234], arr_index[p1*256 + 170], arr_index[p1*256 + 234], p1%2==0);
        compareSwap(arr[p1*256 + 171], arr[p1*256 + 235], arr_index[p1*256 + 171], arr_index[p1*256 + 235], p1%2==0);
        compareSwap(arr[p1*256 + 172], arr[p1*256 + 236], arr_index[p1*256 + 172], arr_index[p1*256 + 236], p1%2==0);
        compareSwap(arr[p1*256 + 173], arr[p1*256 + 237], arr_index[p1*256 + 173], arr_index[p1*256 + 237], p1%2==0);
        compareSwap(arr[p1*256 + 174], arr[p1*256 + 238], arr_index[p1*256 + 174], arr_index[p1*256 + 238], p1%2==0);
        compareSwap(arr[p1*256 + 175], arr[p1*256 + 239], arr_index[p1*256 + 175], arr_index[p1*256 + 239], p1%2==0);
        compareSwap(arr[p1*256 + 176], arr[p1*256 + 240], arr_index[p1*256 + 176], arr_index[p1*256 + 240], p1%2==0);
        compareSwap(arr[p1*256 + 177], arr[p1*256 + 241], arr_index[p1*256 + 177], arr_index[p1*256 + 241], p1%2==0);
        compareSwap(arr[p1*256 + 178], arr[p1*256 + 242], arr_index[p1*256 + 178], arr_index[p1*256 + 242], p1%2==0);
        compareSwap(arr[p1*256 + 179], arr[p1*256 + 243], arr_index[p1*256 + 179], arr_index[p1*256 + 243], p1%2==0);
        compareSwap(arr[p1*256 + 180], arr[p1*256 + 244], arr_index[p1*256 + 180], arr_index[p1*256 + 244], p1%2==0);
        compareSwap(arr[p1*256 + 181], arr[p1*256 + 245], arr_index[p1*256 + 181], arr_index[p1*256 + 245], p1%2==0);
        compareSwap(arr[p1*256 + 182], arr[p1*256 + 246], arr_index[p1*256 + 182], arr_index[p1*256 + 246], p1%2==0);
        compareSwap(arr[p1*256 + 183], arr[p1*256 + 247], arr_index[p1*256 + 183], arr_index[p1*256 + 247], p1%2==0);
        compareSwap(arr[p1*256 + 184], arr[p1*256 + 248], arr_index[p1*256 + 184], arr_index[p1*256 + 248], p1%2==0);
        compareSwap(arr[p1*256 + 185], arr[p1*256 + 249], arr_index[p1*256 + 185], arr_index[p1*256 + 249], p1%2==0);
        compareSwap(arr[p1*256 + 186], arr[p1*256 + 250], arr_index[p1*256 + 186], arr_index[p1*256 + 250], p1%2==0);
        compareSwap(arr[p1*256 + 187], arr[p1*256 + 251], arr_index[p1*256 + 187], arr_index[p1*256 + 251], p1%2==0);
        compareSwap(arr[p1*256 + 188], arr[p1*256 + 252], arr_index[p1*256 + 188], arr_index[p1*256 + 252], p1%2==0);
        compareSwap(arr[p1*256 + 189], arr[p1*256 + 253], arr_index[p1*256 + 189], arr_index[p1*256 + 253], p1%2==0);
        compareSwap(arr[p1*256 + 190], arr[p1*256 + 254], arr_index[p1*256 + 190], arr_index[p1*256 + 254], p1%2==0);
        compareSwap(arr[p1*256 + 191], arr[p1*256 + 255], arr_index[p1*256 + 191], arr_index[p1*256 + 255], p1%2==0);
    }

    for(int p2 = 0; p2 < 1; p2++)
    {
        compareSwap(arr[p2*256 + 0], arr[p2*256 + 32], arr_index[p2*256 + 0], arr_index[p2*256 + 32], p2%2==0);
        compareSwap(arr[p2*256 + 1], arr[p2*256 + 33], arr_index[p2*256 + 1], arr_index[p2*256 + 33], p2%2==0);
        compareSwap(arr[p2*256 + 2], arr[p2*256 + 34], arr_index[p2*256 + 2], arr_index[p2*256 + 34], p2%2==0);
        compareSwap(arr[p2*256 + 3], arr[p2*256 + 35], arr_index[p2*256 + 3], arr_index[p2*256 + 35], p2%2==0);
        compareSwap(arr[p2*256 + 4], arr[p2*256 + 36], arr_index[p2*256 + 4], arr_index[p2*256 + 36], p2%2==0);
        compareSwap(arr[p2*256 + 5], arr[p2*256 + 37], arr_index[p2*256 + 5], arr_index[p2*256 + 37], p2%2==0);
        compareSwap(arr[p2*256 + 6], arr[p2*256 + 38], arr_index[p2*256 + 6], arr_index[p2*256 + 38], p2%2==0);
        compareSwap(arr[p2*256 + 7], arr[p2*256 + 39], arr_index[p2*256 + 7], arr_index[p2*256 + 39], p2%2==0);
        compareSwap(arr[p2*256 + 8], arr[p2*256 + 40], arr_index[p2*256 + 8], arr_index[p2*256 + 40], p2%2==0);
        compareSwap(arr[p2*256 + 9], arr[p2*256 + 41], arr_index[p2*256 + 9], arr_index[p2*256 + 41], p2%2==0);
        compareSwap(arr[p2*256 + 10], arr[p2*256 + 42], arr_index[p2*256 + 10], arr_index[p2*256 + 42], p2%2==0);
        compareSwap(arr[p2*256 + 11], arr[p2*256 + 43], arr_index[p2*256 + 11], arr_index[p2*256 + 43], p2%2==0);
        compareSwap(arr[p2*256 + 12], arr[p2*256 + 44], arr_index[p2*256 + 12], arr_index[p2*256 + 44], p2%2==0);
        compareSwap(arr[p2*256 + 13], arr[p2*256 + 45], arr_index[p2*256 + 13], arr_index[p2*256 + 45], p2%2==0);
        compareSwap(arr[p2*256 + 14], arr[p2*256 + 46], arr_index[p2*256 + 14], arr_index[p2*256 + 46], p2%2==0);
        compareSwap(arr[p2*256 + 15], arr[p2*256 + 47], arr_index[p2*256 + 15], arr_index[p2*256 + 47], p2%2==0);
        compareSwap(arr[p2*256 + 16], arr[p2*256 + 48], arr_index[p2*256 + 16], arr_index[p2*256 + 48], p2%2==0);
        compareSwap(arr[p2*256 + 17], arr[p2*256 + 49], arr_index[p2*256 + 17], arr_index[p2*256 + 49], p2%2==0);
        compareSwap(arr[p2*256 + 18], arr[p2*256 + 50], arr_index[p2*256 + 18], arr_index[p2*256 + 50], p2%2==0);
        compareSwap(arr[p2*256 + 19], arr[p2*256 + 51], arr_index[p2*256 + 19], arr_index[p2*256 + 51], p2%2==0);
        compareSwap(arr[p2*256 + 20], arr[p2*256 + 52], arr_index[p2*256 + 20], arr_index[p2*256 + 52], p2%2==0);
        compareSwap(arr[p2*256 + 21], arr[p2*256 + 53], arr_index[p2*256 + 21], arr_index[p2*256 + 53], p2%2==0);
        compareSwap(arr[p2*256 + 22], arr[p2*256 + 54], arr_index[p2*256 + 22], arr_index[p2*256 + 54], p2%2==0);
        compareSwap(arr[p2*256 + 23], arr[p2*256 + 55], arr_index[p2*256 + 23], arr_index[p2*256 + 55], p2%2==0);
        compareSwap(arr[p2*256 + 24], arr[p2*256 + 56], arr_index[p2*256 + 24], arr_index[p2*256 + 56], p2%2==0);
        compareSwap(arr[p2*256 + 25], arr[p2*256 + 57], arr_index[p2*256 + 25], arr_index[p2*256 + 57], p2%2==0);
        compareSwap(arr[p2*256 + 26], arr[p2*256 + 58], arr_index[p2*256 + 26], arr_index[p2*256 + 58], p2%2==0);
        compareSwap(arr[p2*256 + 27], arr[p2*256 + 59], arr_index[p2*256 + 27], arr_index[p2*256 + 59], p2%2==0);
        compareSwap(arr[p2*256 + 28], arr[p2*256 + 60], arr_index[p2*256 + 28], arr_index[p2*256 + 60], p2%2==0);
        compareSwap(arr[p2*256 + 29], arr[p2*256 + 61], arr_index[p2*256 + 29], arr_index[p2*256 + 61], p2%2==0);
        compareSwap(arr[p2*256 + 30], arr[p2*256 + 62], arr_index[p2*256 + 30], arr_index[p2*256 + 62], p2%2==0);
        compareSwap(arr[p2*256 + 31], arr[p2*256 + 63], arr_index[p2*256 + 31], arr_index[p2*256 + 63], p2%2==0);

        compareSwap(arr[p2*256 + 64], arr[p2*256 + 96], arr_index[p2*256 + 64], arr_index[p2*256 + 96], p2%2==0);
        compareSwap(arr[p2*256 + 65], arr[p2*256 + 97], arr_index[p2*256 + 65], arr_index[p2*256 + 97], p2%2==0);
        compareSwap(arr[p2*256 + 66], arr[p2*256 + 98], arr_index[p2*256 + 66], arr_index[p2*256 + 98], p2%2==0);
        compareSwap(arr[p2*256 + 67], arr[p2*256 + 99], arr_index[p2*256 + 67], arr_index[p2*256 + 99], p2%2==0);
        compareSwap(arr[p2*256 + 68], arr[p2*256 + 100], arr_index[p2*256 + 68], arr_index[p2*256 + 100], p2%2==0);
        compareSwap(arr[p2*256 + 69], arr[p2*256 + 101], arr_index[p2*256 + 69], arr_index[p2*256 + 101], p2%2==0);
        compareSwap(arr[p2*256 + 70], arr[p2*256 + 102], arr_index[p2*256 + 70], arr_index[p2*256 + 102], p2%2==0);
        compareSwap(arr[p2*256 + 71], arr[p2*256 + 103], arr_index[p2*256 + 71], arr_index[p2*256 + 103], p2%2==0);
        compareSwap(arr[p2*256 + 72], arr[p2*256 + 104], arr_index[p2*256 + 72], arr_index[p2*256 + 104], p2%2==0);
        compareSwap(arr[p2*256 + 73], arr[p2*256 + 105], arr_index[p2*256 + 73], arr_index[p2*256 + 105], p2%2==0);
        compareSwap(arr[p2*256 + 74], arr[p2*256 + 106], arr_index[p2*256 + 74], arr_index[p2*256 + 106], p2%2==0);
        compareSwap(arr[p2*256 + 75], arr[p2*256 + 107], arr_index[p2*256 + 75], arr_index[p2*256 + 107], p2%2==0);
        compareSwap(arr[p2*256 + 76], arr[p2*256 + 108], arr_index[p2*256 + 76], arr_index[p2*256 + 108], p2%2==0);
        compareSwap(arr[p2*256 + 77], arr[p2*256 + 109], arr_index[p2*256 + 77], arr_index[p2*256 + 109], p2%2==0);
        compareSwap(arr[p2*256 + 78], arr[p2*256 + 110], arr_index[p2*256 + 78], arr_index[p2*256 + 110], p2%2==0);
        compareSwap(arr[p2*256 + 79], arr[p2*256 + 111], arr_index[p2*256 + 79], arr_index[p2*256 + 111], p2%2==0);
        compareSwap(arr[p2*256 + 80], arr[p2*256 + 112], arr_index[p2*256 + 80], arr_index[p2*256 + 112], p2%2==0);
        compareSwap(arr[p2*256 + 81], arr[p2*256 + 113], arr_index[p2*256 + 81], arr_index[p2*256 + 113], p2%2==0);
        compareSwap(arr[p2*256 + 82], arr[p2*256 + 114], arr_index[p2*256 + 82], arr_index[p2*256 + 114], p2%2==0);
        compareSwap(arr[p2*256 + 83], arr[p2*256 + 115], arr_index[p2*256 + 83], arr_index[p2*256 + 115], p2%2==0);
        compareSwap(arr[p2*256 + 84], arr[p2*256 + 116], arr_index[p2*256 + 84], arr_index[p2*256 + 116], p2%2==0);
        compareSwap(arr[p2*256 + 85], arr[p2*256 + 117], arr_index[p2*256 + 85], arr_index[p2*256 + 117], p2%2==0);
        compareSwap(arr[p2*256 + 86], arr[p2*256 + 118], arr_index[p2*256 + 86], arr_index[p2*256 + 118], p2%2==0);
        compareSwap(arr[p2*256 + 87], arr[p2*256 + 119], arr_index[p2*256 + 87], arr_index[p2*256 + 119], p2%2==0);
        compareSwap(arr[p2*256 + 88], arr[p2*256 + 120], arr_index[p2*256 + 88], arr_index[p2*256 + 120], p2%2==0);
        compareSwap(arr[p2*256 + 89], arr[p2*256 + 121], arr_index[p2*256 + 89], arr_index[p2*256 + 121], p2%2==0);
        compareSwap(arr[p2*256 + 90], arr[p2*256 + 122], arr_index[p2*256 + 90], arr_index[p2*256 + 122], p2%2==0);
        compareSwap(arr[p2*256 + 91], arr[p2*256 + 123], arr_index[p2*256 + 91], arr_index[p2*256 + 123], p2%2==0);
        compareSwap(arr[p2*256 + 92], arr[p2*256 + 124], arr_index[p2*256 + 92], arr_index[p2*256 + 124], p2%2==0);
        compareSwap(arr[p2*256 + 93], arr[p2*256 + 125], arr_index[p2*256 + 93], arr_index[p2*256 + 125], p2%2==0);
        compareSwap(arr[p2*256 + 94], arr[p2*256 + 126], arr_index[p2*256 + 94], arr_index[p2*256 + 126], p2%2==0);
        compareSwap(arr[p2*256 + 95], arr[p2*256 + 127], arr_index[p2*256 + 95], arr_index[p2*256 + 127], p2%2==0);

        compareSwap(arr[p2*256 + 128], arr[p2*256 + 160], arr_index[p2*256 + 128], arr_index[p2*256 + 160], p2%2==0);
        compareSwap(arr[p2*256 + 129], arr[p2*256 + 161], arr_index[p2*256 + 129], arr_index[p2*256 + 161], p2%2==0);
        compareSwap(arr[p2*256 + 130], arr[p2*256 + 162], arr_index[p2*256 + 130], arr_index[p2*256 + 162], p2%2==0);
        compareSwap(arr[p2*256 + 131], arr[p2*256 + 163], arr_index[p2*256 + 131], arr_index[p2*256 + 163], p2%2==0);
        compareSwap(arr[p2*256 + 132], arr[p2*256 + 164], arr_index[p2*256 + 132], arr_index[p2*256 + 164], p2%2==0);
        compareSwap(arr[p2*256 + 133], arr[p2*256 + 165], arr_index[p2*256 + 133], arr_index[p2*256 + 165], p2%2==0);
        compareSwap(arr[p2*256 + 134], arr[p2*256 + 166], arr_index[p2*256 + 134], arr_index[p2*256 + 166], p2%2==0);
        compareSwap(arr[p2*256 + 135], arr[p2*256 + 167], arr_index[p2*256 + 135], arr_index[p2*256 + 167], p2%2==0);
        compareSwap(arr[p2*256 + 136], arr[p2*256 + 168], arr_index[p2*256 + 136], arr_index[p2*256 + 168], p2%2==0);
        compareSwap(arr[p2*256 + 137], arr[p2*256 + 169], arr_index[p2*256 + 137], arr_index[p2*256 + 169], p2%2==0);
        compareSwap(arr[p2*256 + 138], arr[p2*256 + 170], arr_index[p2*256 + 138], arr_index[p2*256 + 170], p2%2==0);
        compareSwap(arr[p2*256 + 139], arr[p2*256 + 171], arr_index[p2*256 + 139], arr_index[p2*256 + 171], p2%2==0);
        compareSwap(arr[p2*256 + 140], arr[p2*256 + 172], arr_index[p2*256 + 140], arr_index[p2*256 + 172], p2%2==0);
        compareSwap(arr[p2*256 + 141], arr[p2*256 + 173], arr_index[p2*256 + 141], arr_index[p2*256 + 173], p2%2==0);
        compareSwap(arr[p2*256 + 142], arr[p2*256 + 174], arr_index[p2*256 + 142], arr_index[p2*256 + 174], p2%2==0);
        compareSwap(arr[p2*256 + 143], arr[p2*256 + 175], arr_index[p2*256 + 143], arr_index[p2*256 + 175], p2%2==0);
        compareSwap(arr[p2*256 + 144], arr[p2*256 + 176], arr_index[p2*256 + 144], arr_index[p2*256 + 176], p2%2==0);
        compareSwap(arr[p2*256 + 145], arr[p2*256 + 177], arr_index[p2*256 + 145], arr_index[p2*256 + 177], p2%2==0);
        compareSwap(arr[p2*256 + 146], arr[p2*256 + 178], arr_index[p2*256 + 146], arr_index[p2*256 + 178], p2%2==0);
        compareSwap(arr[p2*256 + 147], arr[p2*256 + 179], arr_index[p2*256 + 147], arr_index[p2*256 + 179], p2%2==0);
        compareSwap(arr[p2*256 + 148], arr[p2*256 + 180], arr_index[p2*256 + 148], arr_index[p2*256 + 180], p2%2==0);
        compareSwap(arr[p2*256 + 149], arr[p2*256 + 181], arr_index[p2*256 + 149], arr_index[p2*256 + 181], p2%2==0);
        compareSwap(arr[p2*256 + 150], arr[p2*256 + 182], arr_index[p2*256 + 150], arr_index[p2*256 + 182], p2%2==0);
        compareSwap(arr[p2*256 + 151], arr[p2*256 + 183], arr_index[p2*256 + 151], arr_index[p2*256 + 183], p2%2==0);
        compareSwap(arr[p2*256 + 152], arr[p2*256 + 184], arr_index[p2*256 + 152], arr_index[p2*256 + 184], p2%2==0);
        compareSwap(arr[p2*256 + 153], arr[p2*256 + 185], arr_index[p2*256 + 153], arr_index[p2*256 + 185], p2%2==0);
        compareSwap(arr[p2*256 + 154], arr[p2*256 + 186], arr_index[p2*256 + 154], arr_index[p2*256 + 186], p2%2==0);
        compareSwap(arr[p2*256 + 155], arr[p2*256 + 187], arr_index[p2*256 + 155], arr_index[p2*256 + 187], p2%2==0);
        compareSwap(arr[p2*256 + 156], arr[p2*256 + 188], arr_index[p2*256 + 156], arr_index[p2*256 + 188], p2%2==0);
        compareSwap(arr[p2*256 + 157], arr[p2*256 + 189], arr_index[p2*256 + 157], arr_index[p2*256 + 189], p2%2==0);
        compareSwap(arr[p2*256 + 158], arr[p2*256 + 190], arr_index[p2*256 + 158], arr_index[p2*256 + 190], p2%2==0);
        compareSwap(arr[p2*256 + 159], arr[p2*256 + 191], arr_index[p2*256 + 159], arr_index[p2*256 + 191], p2%2==0);
        
        compareSwap(arr[p2*256 + 192], arr[p2*256 + 224], arr_index[p2*256 + 192], arr_index[p2*256 + 224], p2%2==0);
        compareSwap(arr[p2*256 + 193], arr[p2*256 + 225], arr_index[p2*256 + 193], arr_index[p2*256 + 225], p2%2==0);
        compareSwap(arr[p2*256 + 194], arr[p2*256 + 226], arr_index[p2*256 + 194], arr_index[p2*256 + 226], p2%2==0);
        compareSwap(arr[p2*256 + 195], arr[p2*256 + 227], arr_index[p2*256 + 195], arr_index[p2*256 + 227], p2%2==0);
        compareSwap(arr[p2*256 + 196], arr[p2*256 + 228], arr_index[p2*256 + 196], arr_index[p2*256 + 228], p2%2==0);
        compareSwap(arr[p2*256 + 197], arr[p2*256 + 229], arr_index[p2*256 + 197], arr_index[p2*256 + 229], p2%2==0);
        compareSwap(arr[p2*256 + 198], arr[p2*256 + 230], arr_index[p2*256 + 198], arr_index[p2*256 + 230], p2%2==0);
        compareSwap(arr[p2*256 + 199], arr[p2*256 + 231], arr_index[p2*256 + 199], arr_index[p2*256 + 231], p2%2==0);
        compareSwap(arr[p2*256 + 200], arr[p2*256 + 232], arr_index[p2*256 + 200], arr_index[p2*256 + 232], p2%2==0);
        compareSwap(arr[p2*256 + 201], arr[p2*256 + 233], arr_index[p2*256 + 201], arr_index[p2*256 + 233], p2%2==0);
        compareSwap(arr[p2*256 + 202], arr[p2*256 + 234], arr_index[p2*256 + 202], arr_index[p2*256 + 234], p2%2==0);
        compareSwap(arr[p2*256 + 203], arr[p2*256 + 235], arr_index[p2*256 + 203], arr_index[p2*256 + 235], p2%2==0);
        compareSwap(arr[p2*256 + 204], arr[p2*256 + 236], arr_index[p2*256 + 204], arr_index[p2*256 + 236], p2%2==0);
        compareSwap(arr[p2*256 + 205], arr[p2*256 + 237], arr_index[p2*256 + 205], arr_index[p2*256 + 237], p2%2==0);
        compareSwap(arr[p2*256 + 206], arr[p2*256 + 238], arr_index[p2*256 + 206], arr_index[p2*256 + 238], p2%2==0);
        compareSwap(arr[p2*256 + 207], arr[p2*256 + 239], arr_index[p2*256 + 207], arr_index[p2*256 + 239], p2%2==0);
        compareSwap(arr[p2*256 + 208], arr[p2*256 + 240], arr_index[p2*256 + 208], arr_index[p2*256 + 240], p2%2==0);
        compareSwap(arr[p2*256 + 209], arr[p2*256 + 241], arr_index[p2*256 + 209], arr_index[p2*256 + 241], p2%2==0);
        compareSwap(arr[p2*256 + 210], arr[p2*256 + 242], arr_index[p2*256 + 210], arr_index[p2*256 + 242], p2%2==0);
        compareSwap(arr[p2*256 + 211], arr[p2*256 + 243], arr_index[p2*256 + 211], arr_index[p2*256 + 243], p2%2==0);
        compareSwap(arr[p2*256 + 212], arr[p2*256 + 244], arr_index[p2*256 + 212], arr_index[p2*256 + 244], p2%2==0);
        compareSwap(arr[p2*256 + 213], arr[p2*256 + 245], arr_index[p2*256 + 213], arr_index[p2*256 + 245], p2%2==0);
        compareSwap(arr[p2*256 + 214], arr[p2*256 + 246], arr_index[p2*256 + 214], arr_index[p2*256 + 246], p2%2==0);
        compareSwap(arr[p2*256 + 215], arr[p2*256 + 247], arr_index[p2*256 + 215], arr_index[p2*256 + 247], p2%2==0);
        compareSwap(arr[p2*256 + 216], arr[p2*256 + 248], arr_index[p2*256 + 216], arr_index[p2*256 + 248], p2%2==0);
        compareSwap(arr[p2*256 + 217], arr[p2*256 + 249], arr_index[p2*256 + 217], arr_index[p2*256 + 249], p2%2==0);
        compareSwap(arr[p2*256 + 218], arr[p2*256 + 250], arr_index[p2*256 + 218], arr_index[p2*256 + 250], p2%2==0);
        compareSwap(arr[p2*256 + 219], arr[p2*256 + 251], arr_index[p2*256 + 219], arr_index[p2*256 + 251], p2%2==0);
        compareSwap(arr[p2*256 + 220], arr[p2*256 + 252], arr_index[p2*256 + 220], arr_index[p2*256 + 252], p2%2==0);
        compareSwap(arr[p2*256 + 221], arr[p2*256 + 253], arr_index[p2*256 + 221], arr_index[p2*256 + 253], p2%2==0);
        compareSwap(arr[p2*256 + 222], arr[p2*256 + 254], arr_index[p2*256 + 222], arr_index[p2*256 + 254], p2%2==0);
        compareSwap(arr[p2*256 + 223], arr[p2*256 + 255], arr_index[p2*256 + 223], arr_index[p2*256 + 255], p2%2==0);
    }

    for(int p3 = 0; p3 < 1; p3++)
    {
        compareSwap(arr[p3*256 + 0], arr[p3*256 + 16], arr_index[p3*256 + 0], arr_index[p3*256 + 16], p3%2==0);
        compareSwap(arr[p3*256 + 1], arr[p3*256 + 17], arr_index[p3*256 + 1], arr_index[p3*256 + 17], p3%2==0);
        compareSwap(arr[p3*256 + 2], arr[p3*256 + 18], arr_index[p3*256 + 2], arr_index[p3*256 + 18], p3%2==0);
        compareSwap(arr[p3*256 + 3], arr[p3*256 + 19], arr_index[p3*256 + 3], arr_index[p3*256 + 19], p3%2==0);
        compareSwap(arr[p3*256 + 4], arr[p3*256 + 20], arr_index[p3*256 + 4], arr_index[p3*256 + 20], p3%2==0);
        compareSwap(arr[p3*256 + 5], arr[p3*256 + 21], arr_index[p3*256 + 5], arr_index[p3*256 + 21], p3%2==0);
        compareSwap(arr[p3*256 + 6], arr[p3*256 + 22], arr_index[p3*256 + 6], arr_index[p3*256 + 22], p3%2==0);
        compareSwap(arr[p3*256 + 7], arr[p3*256 + 23], arr_index[p3*256 + 7], arr_index[p3*256 + 23], p3%2==0);
        compareSwap(arr[p3*256 + 8], arr[p3*256 + 24], arr_index[p3*256 + 8], arr_index[p3*256 + 24], p3%2==0);
        compareSwap(arr[p3*256 + 9], arr[p3*256 + 25], arr_index[p3*256 + 9], arr_index[p3*256 + 25], p3%2==0);
        compareSwap(arr[p3*256 + 10], arr[p3*256 + 26], arr_index[p3*256 + 10], arr_index[p3*256 + 26], p3%2==0);
        compareSwap(arr[p3*256 + 11], arr[p3*256 + 27], arr_index[p3*256 + 11], arr_index[p3*256 + 27], p3%2==0);
        compareSwap(arr[p3*256 + 12], arr[p3*256 + 28], arr_index[p3*256 + 12], arr_index[p3*256 + 28], p3%2==0);
        compareSwap(arr[p3*256 + 13], arr[p3*256 + 29], arr_index[p3*256 + 13], arr_index[p3*256 + 29], p3%2==0);
        compareSwap(arr[p3*256 + 14], arr[p3*256 + 30], arr_index[p3*256 + 14], arr_index[p3*256 + 30], p3%2==0);
        compareSwap(arr[p3*256 + 15], arr[p3*256 + 31], arr_index[p3*256 + 15], arr_index[p3*256 + 31], p3%2==0);

        compareSwap(arr[p3*256 + 32], arr[p3*256 + 48], arr_index[p3*256 + 32], arr_index[p3*256 + 48], p3%2==0);
        compareSwap(arr[p3*256 + 33], arr[p3*256 + 49], arr_index[p3*256 + 33], arr_index[p3*256 + 49], p3%2==0);
        compareSwap(arr[p3*256 + 34], arr[p3*256 + 50], arr_index[p3*256 + 34], arr_index[p3*256 + 50], p3%2==0);
        compareSwap(arr[p3*256 + 35], arr[p3*256 + 51], arr_index[p3*256 + 35], arr_index[p3*256 + 51], p3%2==0);
        compareSwap(arr[p3*256 + 36], arr[p3*256 + 52], arr_index[p3*256 + 36], arr_index[p3*256 + 52], p3%2==0);
        compareSwap(arr[p3*256 + 37], arr[p3*256 + 53], arr_index[p3*256 + 37], arr_index[p3*256 + 53], p3%2==0);
        compareSwap(arr[p3*256 + 38], arr[p3*256 + 54], arr_index[p3*256 + 38], arr_index[p3*256 + 54], p3%2==0);
        compareSwap(arr[p3*256 + 39], arr[p3*256 + 55], arr_index[p3*256 + 39], arr_index[p3*256 + 55], p3%2==0);
        compareSwap(arr[p3*256 + 40], arr[p3*256 + 56], arr_index[p3*256 + 40], arr_index[p3*256 + 56], p3%2==0);
        compareSwap(arr[p3*256 + 41], arr[p3*256 + 57], arr_index[p3*256 + 41], arr_index[p3*256 + 57], p3%2==0);
        compareSwap(arr[p3*256 + 42], arr[p3*256 + 58], arr_index[p3*256 + 42], arr_index[p3*256 + 58], p3%2==0);
        compareSwap(arr[p3*256 + 43], arr[p3*256 + 59], arr_index[p3*256 + 43], arr_index[p3*256 + 59], p3%2==0);
        compareSwap(arr[p3*256 + 44], arr[p3*256 + 60], arr_index[p3*256 + 44], arr_index[p3*256 + 60], p3%2==0);
        compareSwap(arr[p3*256 + 45], arr[p3*256 + 61], arr_index[p3*256 + 45], arr_index[p3*256 + 61], p3%2==0);
        compareSwap(arr[p3*256 + 46], arr[p3*256 + 62], arr_index[p3*256 + 46], arr_index[p3*256 + 62], p3%2==0);
        compareSwap(arr[p3*256 + 47], arr[p3*256 + 63], arr_index[p3*256 + 47], arr_index[p3*256 + 63], p3%2==0);

        compareSwap(arr[p3*256 + 64], arr[p3*256 + 80], arr_index[p3*256 + 64], arr_index[p3*256 + 80], p3%2==0);
        compareSwap(arr[p3*256 + 65], arr[p3*256 + 81], arr_index[p3*256 + 65], arr_index[p3*256 + 81], p3%2==0);
        compareSwap(arr[p3*256 + 66], arr[p3*256 + 82], arr_index[p3*256 + 66], arr_index[p3*256 + 82], p3%2==0);
        compareSwap(arr[p3*256 + 67], arr[p3*256 + 83], arr_index[p3*256 + 67], arr_index[p3*256 + 83], p3%2==0);
        compareSwap(arr[p3*256 + 68], arr[p3*256 + 84], arr_index[p3*256 + 68], arr_index[p3*256 + 84], p3%2==0);
        compareSwap(arr[p3*256 + 69], arr[p3*256 + 85], arr_index[p3*256 + 69], arr_index[p3*256 + 85], p3%2==0);
        compareSwap(arr[p3*256 + 70], arr[p3*256 + 86], arr_index[p3*256 + 70], arr_index[p3*256 + 86], p3%2==0);
        compareSwap(arr[p3*256 + 71], arr[p3*256 + 87], arr_index[p3*256 + 71], arr_index[p3*256 + 87], p3%2==0);
        compareSwap(arr[p3*256 + 72], arr[p3*256 + 88], arr_index[p3*256 + 72], arr_index[p3*256 + 88], p3%2==0);
        compareSwap(arr[p3*256 + 73], arr[p3*256 + 89], arr_index[p3*256 + 73], arr_index[p3*256 + 89], p3%2==0);
        compareSwap(arr[p3*256 + 74], arr[p3*256 + 90], arr_index[p3*256 + 74], arr_index[p3*256 + 90], p3%2==0);
        compareSwap(arr[p3*256 + 75], arr[p3*256 + 91], arr_index[p3*256 + 75], arr_index[p3*256 + 91], p3%2==0);
        compareSwap(arr[p3*256 + 76], arr[p3*256 + 92], arr_index[p3*256 + 76], arr_index[p3*256 + 92], p3%2==0);
        compareSwap(arr[p3*256 + 77], arr[p3*256 + 93], arr_index[p3*256 + 77], arr_index[p3*256 + 93], p3%2==0);
        compareSwap(arr[p3*256 + 78], arr[p3*256 + 94], arr_index[p3*256 + 78], arr_index[p3*256 + 94], p3%2==0);
        compareSwap(arr[p3*256 + 79], arr[p3*256 + 95], arr_index[p3*256 + 79], arr_index[p3*256 + 95], p3%2==0);

        compareSwap(arr[p3*256 + 96], arr[p3*256 + 112], arr_index[p3*256 + 96], arr_index[p3*256 + 112], p3%2==0);
        compareSwap(arr[p3*256 + 97], arr[p3*256 + 113], arr_index[p3*256 + 97], arr_index[p3*256 + 113], p3%2==0);
        compareSwap(arr[p3*256 + 98], arr[p3*256 + 114], arr_index[p3*256 + 98], arr_index[p3*256 + 114], p3%2==0);
        compareSwap(arr[p3*256 + 99], arr[p3*256 + 115], arr_index[p3*256 + 99], arr_index[p3*256 + 115], p3%2==0);
        compareSwap(arr[p3*256 + 100], arr[p3*256 + 116], arr_index[p3*256 + 100], arr_index[p3*256 + 116], p3%2==0);
        compareSwap(arr[p3*256 + 101], arr[p3*256 + 117], arr_index[p3*256 + 101], arr_index[p3*256 + 117], p3%2==0);
        compareSwap(arr[p3*256 + 102], arr[p3*256 + 118], arr_index[p3*256 + 102], arr_index[p3*256 + 118], p3%2==0);
        compareSwap(arr[p3*256 + 103], arr[p3*256 + 119], arr_index[p3*256 + 103], arr_index[p3*256 + 119], p3%2==0);
        compareSwap(arr[p3*256 + 104], arr[p3*256 + 120], arr_index[p3*256 + 104], arr_index[p3*256 + 120], p3%2==0);
        compareSwap(arr[p3*256 + 105], arr[p3*256 + 121], arr_index[p3*256 + 105], arr_index[p3*256 + 121], p3%2==0);
        compareSwap(arr[p3*256 + 106], arr[p3*256 + 122], arr_index[p3*256 + 106], arr_index[p3*256 + 122], p3%2==0);
        compareSwap(arr[p3*256 + 107], arr[p3*256 + 123], arr_index[p3*256 + 107], arr_index[p3*256 + 123], p3%2==0);
        compareSwap(arr[p3*256 + 108], arr[p3*256 + 124], arr_index[p3*256 + 108], arr_index[p3*256 + 124], p3%2==0);
        compareSwap(arr[p3*256 + 109], arr[p3*256 + 125], arr_index[p3*256 + 109], arr_index[p3*256 + 125], p3%2==0);
        compareSwap(arr[p3*256 + 110], arr[p3*256 + 126], arr_index[p3*256 + 110], arr_index[p3*256 + 126], p3%2==0);
        compareSwap(arr[p3*256 + 111], arr[p3*256 + 127], arr_index[p3*256 + 111], arr_index[p3*256 + 127], p3%2==0);

        compareSwap(arr[p3*256 + 128], arr[p3*256 + 144], arr_index[p3*256 + 128], arr_index[p3*256 + 144], p3%2==0);
        compareSwap(arr[p3*256 + 129], arr[p3*256 + 145], arr_index[p3*256 + 129], arr_index[p3*256 + 145], p3%2==0);
        compareSwap(arr[p3*256 + 130], arr[p3*256 + 146], arr_index[p3*256 + 130], arr_index[p3*256 + 146], p3%2==0);
        compareSwap(arr[p3*256 + 131], arr[p3*256 + 147], arr_index[p3*256 + 131], arr_index[p3*256 + 147], p3%2==0);
        compareSwap(arr[p3*256 + 132], arr[p3*256 + 148], arr_index[p3*256 + 132], arr_index[p3*256 + 148], p3%2==0);
        compareSwap(arr[p3*256 + 133], arr[p3*256 + 149], arr_index[p3*256 + 133], arr_index[p3*256 + 149], p3%2==0);
        compareSwap(arr[p3*256 + 134], arr[p3*256 + 150], arr_index[p3*256 + 134], arr_index[p3*256 + 150], p3%2==0);
        compareSwap(arr[p3*256 + 135], arr[p3*256 + 151], arr_index[p3*256 + 135], arr_index[p3*256 + 151], p3%2==0);
        compareSwap(arr[p3*256 + 136], arr[p3*256 + 152], arr_index[p3*256 + 136], arr_index[p3*256 + 152], p3%2==0);
        compareSwap(arr[p3*256 + 137], arr[p3*256 + 153], arr_index[p3*256 + 137], arr_index[p3*256 + 153], p3%2==0);
        compareSwap(arr[p3*256 + 138], arr[p3*256 + 154], arr_index[p3*256 + 138], arr_index[p3*256 + 154], p3%2==0);
        compareSwap(arr[p3*256 + 139], arr[p3*256 + 155], arr_index[p3*256 + 139], arr_index[p3*256 + 155], p3%2==0);
        compareSwap(arr[p3*256 + 140], arr[p3*256 + 156], arr_index[p3*256 + 140], arr_index[p3*256 + 156], p3%2==0);
        compareSwap(arr[p3*256 + 141], arr[p3*256 + 157], arr_index[p3*256 + 141], arr_index[p3*256 + 157], p3%2==0);
        compareSwap(arr[p3*256 + 142], arr[p3*256 + 158], arr_index[p3*256 + 142], arr_index[p3*256 + 158], p3%2==0);
        compareSwap(arr[p3*256 + 143], arr[p3*256 + 159], arr_index[p3*256 + 143], arr_index[p3*256 + 159], p3%2==0);

        compareSwap(arr[p3*256 + 160], arr[p3*256 + 176], arr_index[p3*256 + 160], arr_index[p3*256 + 176], p3%2==0);
        compareSwap(arr[p3*256 + 161], arr[p3*256 + 177], arr_index[p3*256 + 161], arr_index[p3*256 + 177], p3%2==0);
        compareSwap(arr[p3*256 + 162], arr[p3*256 + 178], arr_index[p3*256 + 162], arr_index[p3*256 + 178], p3%2==0);
        compareSwap(arr[p3*256 + 163], arr[p3*256 + 179], arr_index[p3*256 + 163], arr_index[p3*256 + 179], p3%2==0);
        compareSwap(arr[p3*256 + 164], arr[p3*256 + 180], arr_index[p3*256 + 164], arr_index[p3*256 + 180], p3%2==0);
        compareSwap(arr[p3*256 + 165], arr[p3*256 + 181], arr_index[p3*256 + 165], arr_index[p3*256 + 181], p3%2==0);
        compareSwap(arr[p3*256 + 166], arr[p3*256 + 182], arr_index[p3*256 + 166], arr_index[p3*256 + 182], p3%2==0);
        compareSwap(arr[p3*256 + 167], arr[p3*256 + 183], arr_index[p3*256 + 167], arr_index[p3*256 + 183], p3%2==0);
        compareSwap(arr[p3*256 + 168], arr[p3*256 + 184], arr_index[p3*256 + 168], arr_index[p3*256 + 184], p3%2==0);
        compareSwap(arr[p3*256 + 169], arr[p3*256 + 185], arr_index[p3*256 + 169], arr_index[p3*256 + 185], p3%2==0);
        compareSwap(arr[p3*256 + 170], arr[p3*256 + 186], arr_index[p3*256 + 170], arr_index[p3*256 + 186], p3%2==0);
        compareSwap(arr[p3*256 + 171], arr[p3*256 + 187], arr_index[p3*256 + 171], arr_index[p3*256 + 187], p3%2==0);
        compareSwap(arr[p3*256 + 172], arr[p3*256 + 188], arr_index[p3*256 + 172], arr_index[p3*256 + 188], p3%2==0);
        compareSwap(arr[p3*256 + 173], arr[p3*256 + 189], arr_index[p3*256 + 173], arr_index[p3*256 + 189], p3%2==0);
        compareSwap(arr[p3*256 + 174], arr[p3*256 + 190], arr_index[p3*256 + 174], arr_index[p3*256 + 190], p3%2==0);
        compareSwap(arr[p3*256 + 175], arr[p3*256 + 191], arr_index[p3*256 + 175], arr_index[p3*256 + 191], p3%2==0);

        compareSwap(arr[p3*256 + 192], arr[p3*256 + 208], arr_index[p3*256 + 192], arr_index[p3*256 + 208], p3%2==0);
        compareSwap(arr[p3*256 + 193], arr[p3*256 + 209], arr_index[p3*256 + 193], arr_index[p3*256 + 209], p3%2==0);
        compareSwap(arr[p3*256 + 194], arr[p3*256 + 210], arr_index[p3*256 + 194], arr_index[p3*256 + 210], p3%2==0);
        compareSwap(arr[p3*256 + 195], arr[p3*256 + 211], arr_index[p3*256 + 195], arr_index[p3*256 + 211], p3%2==0);
        compareSwap(arr[p3*256 + 196], arr[p3*256 + 212], arr_index[p3*256 + 196], arr_index[p3*256 + 212], p3%2==0);
        compareSwap(arr[p3*256 + 197], arr[p3*256 + 213], arr_index[p3*256 + 197], arr_index[p3*256 + 213], p3%2==0);
        compareSwap(arr[p3*256 + 198], arr[p3*256 + 214], arr_index[p3*256 + 198], arr_index[p3*256 + 214], p3%2==0);
        compareSwap(arr[p3*256 + 199], arr[p3*256 + 215], arr_index[p3*256 + 199], arr_index[p3*256 + 215], p3%2==0);
        compareSwap(arr[p3*256 + 200], arr[p3*256 + 216], arr_index[p3*256 + 200], arr_index[p3*256 + 216], p3%2==0);
        compareSwap(arr[p3*256 + 201], arr[p3*256 + 217], arr_index[p3*256 + 201], arr_index[p3*256 + 217], p3%2==0);
        compareSwap(arr[p3*256 + 202], arr[p3*256 + 218], arr_index[p3*256 + 202], arr_index[p3*256 + 218], p3%2==0);
        compareSwap(arr[p3*256 + 203], arr[p3*256 + 219], arr_index[p3*256 + 203], arr_index[p3*256 + 219], p3%2==0);
        compareSwap(arr[p3*256 + 204], arr[p3*256 + 220], arr_index[p3*256 + 204], arr_index[p3*256 + 220], p3%2==0);
        compareSwap(arr[p3*256 + 205], arr[p3*256 + 221], arr_index[p3*256 + 205], arr_index[p3*256 + 221], p3%2==0);
        compareSwap(arr[p3*256 + 206], arr[p3*256 + 222], arr_index[p3*256 + 206], arr_index[p3*256 + 222], p3%2==0);
        compareSwap(arr[p3*256 + 207], arr[p3*256 + 223], arr_index[p3*256 + 207], arr_index[p3*256 + 223], p3%2==0);

        compareSwap(arr[p3*256 + 224], arr[p3*256 + 240], arr_index[p3*256 + 224], arr_index[p3*256 + 240], p3%2==0);
        compareSwap(arr[p3*256 + 225], arr[p3*256 + 241], arr_index[p3*256 + 225], arr_index[p3*256 + 241], p3%2==0);
        compareSwap(arr[p3*256 + 226], arr[p3*256 + 242], arr_index[p3*256 + 226], arr_index[p3*256 + 242], p3%2==0);
        compareSwap(arr[p3*256 + 227], arr[p3*256 + 243], arr_index[p3*256 + 227], arr_index[p3*256 + 243], p3%2==0);
        compareSwap(arr[p3*256 + 228], arr[p3*256 + 244], arr_index[p3*256 + 228], arr_index[p3*256 + 244], p3%2==0);
        compareSwap(arr[p3*256 + 229], arr[p3*256 + 245], arr_index[p3*256 + 229], arr_index[p3*256 + 245], p3%2==0);
        compareSwap(arr[p3*256 + 230], arr[p3*256 + 246], arr_index[p3*256 + 230], arr_index[p3*256 + 246], p3%2==0);
        compareSwap(arr[p3*256 + 231], arr[p3*256 + 247], arr_index[p3*256 + 231], arr_index[p3*256 + 247], p3%2==0);
        compareSwap(arr[p3*256 + 232], arr[p3*256 + 248], arr_index[p3*256 + 232], arr_index[p3*256 + 248], p3%2==0);
        compareSwap(arr[p3*256 + 233], arr[p3*256 + 249], arr_index[p3*256 + 233], arr_index[p3*256 + 249], p3%2==0);
        compareSwap(arr[p3*256 + 234], arr[p3*256 + 250], arr_index[p3*256 + 234], arr_index[p3*256 + 250], p3%2==0);
        compareSwap(arr[p3*256 + 235], arr[p3*256 + 251], arr_index[p3*256 + 235], arr_index[p3*256 + 251], p3%2==0);
        compareSwap(arr[p3*256 + 236], arr[p3*256 + 252], arr_index[p3*256 + 236], arr_index[p3*256 + 252], p3%2==0);
        compareSwap(arr[p3*256 + 237], arr[p3*256 + 253], arr_index[p3*256 + 237], arr_index[p3*256 + 253], p3%2==0);
        compareSwap(arr[p3*256 + 238], arr[p3*256 + 254], arr_index[p3*256 + 238], arr_index[p3*256 + 254], p3%2==0);
        compareSwap(arr[p3*256 + 239], arr[p3*256 + 255], arr_index[p3*256 + 239], arr_index[p3*256 + 255], p3%2==0);
    }

    for(int p4 = 0; p4 < 1; p4++)
    {
        compareSwap(arr[p4*256 + 0], arr[p4*256 + 8], arr_index[p4*256 + 0], arr_index[p4*256 + 8], p4%2==0);
        compareSwap(arr[p4*256 + 1], arr[p4*256 + 9], arr_index[p4*256 + 1], arr_index[p4*256 + 9], p4%2==0);
        compareSwap(arr[p4*256 + 2], arr[p4*256 + 10], arr_index[p4*256 + 2], arr_index[p4*256 + 10], p4%2==0);
        compareSwap(arr[p4*256 + 3], arr[p4*256 + 11], arr_index[p4*256 + 3], arr_index[p4*256 + 11], p4%2==0);
        compareSwap(arr[p4*256 + 4], arr[p4*256 + 12], arr_index[p4*256 + 4], arr_index[p4*256 + 12], p4%2==0);
        compareSwap(arr[p4*256 + 5], arr[p4*256 + 13], arr_index[p4*256 + 5], arr_index[p4*256 + 13], p4%2==0);
        compareSwap(arr[p4*256 + 6], arr[p4*256 + 14], arr_index[p4*256 + 6], arr_index[p4*256 + 14], p4%2==0);
        compareSwap(arr[p4*256 + 7], arr[p4*256 + 15], arr_index[p4*256 + 7], arr_index[p4*256 + 15], p4%2==0);

        compareSwap(arr[p4*256 + 16], arr[p4*256 + 24], arr_index[p4*256 + 16], arr_index[p4*256 + 24], p4%2==0);
        compareSwap(arr[p4*256 + 17], arr[p4*256 + 25], arr_index[p4*256 + 17], arr_index[p4*256 + 25], p4%2==0);
        compareSwap(arr[p4*256 + 18], arr[p4*256 + 26], arr_index[p4*256 + 18], arr_index[p4*256 + 26], p4%2==0);
        compareSwap(arr[p4*256 + 19], arr[p4*256 + 27], arr_index[p4*256 + 19], arr_index[p4*256 + 27], p4%2==0);
        compareSwap(arr[p4*256 + 20], arr[p4*256 + 28], arr_index[p4*256 + 20], arr_index[p4*256 + 28], p4%2==0);
        compareSwap(arr[p4*256 + 21], arr[p4*256 + 29], arr_index[p4*256 + 21], arr_index[p4*256 + 29], p4%2==0);
        compareSwap(arr[p4*256 + 22], arr[p4*256 + 30], arr_index[p4*256 + 22], arr_index[p4*256 + 30], p4%2==0);
        compareSwap(arr[p4*256 + 23], arr[p4*256 + 31], arr_index[p4*256 + 23], arr_index[p4*256 + 31], p4%2==0);

        compareSwap(arr[p4*256 + 32], arr[p4*256 + 40], arr_index[p4*256 + 32], arr_index[p4*256 + 40], p4%2==0);
        compareSwap(arr[p4*256 + 33], arr[p4*256 + 41], arr_index[p4*256 + 33], arr_index[p4*256 + 41], p4%2==0);
        compareSwap(arr[p4*256 + 34], arr[p4*256 + 42], arr_index[p4*256 + 34], arr_index[p4*256 + 42], p4%2==0);
        compareSwap(arr[p4*256 + 35], arr[p4*256 + 43], arr_index[p4*256 + 35], arr_index[p4*256 + 43], p4%2==0);
        compareSwap(arr[p4*256 + 36], arr[p4*256 + 44], arr_index[p4*256 + 36], arr_index[p4*256 + 44], p4%2==0);
        compareSwap(arr[p4*256 + 37], arr[p4*256 + 45], arr_index[p4*256 + 37], arr_index[p4*256 + 45], p4%2==0);
        compareSwap(arr[p4*256 + 38], arr[p4*256 + 46], arr_index[p4*256 + 38], arr_index[p4*256 + 46], p4%2==0);
        compareSwap(arr[p4*256 + 39], arr[p4*256 + 47], arr_index[p4*256 + 39], arr_index[p4*256 + 47], p4%2==0);

        compareSwap(arr[p4*256 + 48], arr[p4*256 + 56], arr_index[p4*256 + 48], arr_index[p4*256 + 56], p4%2==0);
        compareSwap(arr[p4*256 + 49], arr[p4*256 + 57], arr_index[p4*256 + 49], arr_index[p4*256 + 57], p4%2==0);
        compareSwap(arr[p4*256 + 50], arr[p4*256 + 58], arr_index[p4*256 + 50], arr_index[p4*256 + 58], p4%2==0);
        compareSwap(arr[p4*256 + 51], arr[p4*256 + 59], arr_index[p4*256 + 51], arr_index[p4*256 + 59], p4%2==0);
        compareSwap(arr[p4*256 + 52], arr[p4*256 + 60], arr_index[p4*256 + 52], arr_index[p4*256 + 60], p4%2==0);
        compareSwap(arr[p4*256 + 53], arr[p4*256 + 61], arr_index[p4*256 + 53], arr_index[p4*256 + 61], p4%2==0);
        compareSwap(arr[p4*256 + 54], arr[p4*256 + 62], arr_index[p4*256 + 54], arr_index[p4*256 + 62], p4%2==0);
        compareSwap(arr[p4*256 + 55], arr[p4*256 + 63], arr_index[p4*256 + 55], arr_index[p4*256 + 63], p4%2==0);

        compareSwap(arr[p4*256 + 64], arr[p4*256 + 72], arr_index[p4*256 + 64], arr_index[p4*256 + 72], p4%2==0);
        compareSwap(arr[p4*256 + 65], arr[p4*256 + 73], arr_index[p4*256 + 65], arr_index[p4*256 + 73], p4%2==0);
        compareSwap(arr[p4*256 + 66], arr[p4*256 + 74], arr_index[p4*256 + 66], arr_index[p4*256 + 74], p4%2==0);
        compareSwap(arr[p4*256 + 67], arr[p4*256 + 75], arr_index[p4*256 + 67], arr_index[p4*256 + 75], p4%2==0);
        compareSwap(arr[p4*256 + 68], arr[p4*256 + 76], arr_index[p4*256 + 68], arr_index[p4*256 + 76], p4%2==0);
        compareSwap(arr[p4*256 + 69], arr[p4*256 + 77], arr_index[p4*256 + 69], arr_index[p4*256 + 77], p4%2==0);
        compareSwap(arr[p4*256 + 70], arr[p4*256 + 78], arr_index[p4*256 + 70], arr_index[p4*256 + 78], p4%2==0);
        compareSwap(arr[p4*256 + 71], arr[p4*256 + 79], arr_index[p4*256 + 71], arr_index[p4*256 + 79], p4%2==0);

        compareSwap(arr[p4*256 + 80], arr[p4*256 + 88], arr_index[p4*256 + 80], arr_index[p4*256 + 88], p4%2==0);
        compareSwap(arr[p4*256 + 81], arr[p4*256 + 89], arr_index[p4*256 + 81], arr_index[p4*256 + 89], p4%2==0);
        compareSwap(arr[p4*256 + 82], arr[p4*256 + 90], arr_index[p4*256 + 82], arr_index[p4*256 + 90], p4%2==0);
        compareSwap(arr[p4*256 + 83], arr[p4*256 + 91], arr_index[p4*256 + 83], arr_index[p4*256 + 91], p4%2==0);
        compareSwap(arr[p4*256 + 84], arr[p4*256 + 92], arr_index[p4*256 + 84], arr_index[p4*256 + 92], p4%2==0);
        compareSwap(arr[p4*256 + 85], arr[p4*256 + 93], arr_index[p4*256 + 85], arr_index[p4*256 + 93], p4%2==0);
        compareSwap(arr[p4*256 + 86], arr[p4*256 + 94], arr_index[p4*256 + 86], arr_index[p4*256 + 94], p4%2==0);
        compareSwap(arr[p4*256 + 87], arr[p4*256 + 95], arr_index[p4*256 + 87], arr_index[p4*256 + 95], p4%2==0);

        compareSwap(arr[p4*256 + 96], arr[p4*256 + 104], arr_index[p4*256 + 96], arr_index[p4*256 + 104], p4%2==0);
        compareSwap(arr[p4*256 + 97], arr[p4*256 + 105], arr_index[p4*256 + 97], arr_index[p4*256 + 105], p4%2==0);
        compareSwap(arr[p4*256 + 98], arr[p4*256 + 106], arr_index[p4*256 + 98], arr_index[p4*256 + 106], p4%2==0);
        compareSwap(arr[p4*256 + 99], arr[p4*256 + 107], arr_index[p4*256 + 99], arr_index[p4*256 + 107], p4%2==0);
        compareSwap(arr[p4*256 + 100], arr[p4*256 + 108], arr_index[p4*256 + 100], arr_index[p4*256 + 108], p4%2==0);
        compareSwap(arr[p4*256 + 101], arr[p4*256 + 109], arr_index[p4*256 + 101], arr_index[p4*256 + 109], p4%2==0);
        compareSwap(arr[p4*256 + 102], arr[p4*256 + 110], arr_index[p4*256 + 102], arr_index[p4*256 + 110], p4%2==0);
        compareSwap(arr[p4*256 + 103], arr[p4*256 + 111], arr_index[p4*256 + 103], arr_index[p4*256 + 111], p4%2==0);

        compareSwap(arr[p4*256 + 112], arr[p4*256 + 120], arr_index[p4*256 + 112], arr_index[p4*256 + 120], p4%2==0);
        compareSwap(arr[p4*256 + 113], arr[p4*256 + 121], arr_index[p4*256 + 113], arr_index[p4*256 + 121], p4%2==0);
        compareSwap(arr[p4*256 + 114], arr[p4*256 + 122], arr_index[p4*256 + 114], arr_index[p4*256 + 122], p4%2==0);
        compareSwap(arr[p4*256 + 115], arr[p4*256 + 123], arr_index[p4*256 + 115], arr_index[p4*256 + 123], p4%2==0);
        compareSwap(arr[p4*256 + 116], arr[p4*256 + 124], arr_index[p4*256 + 116], arr_index[p4*256 + 124], p4%2==0);
        compareSwap(arr[p4*256 + 117], arr[p4*256 + 125], arr_index[p4*256 + 117], arr_index[p4*256 + 125], p4%2==0);
        compareSwap(arr[p4*256 + 118], arr[p4*256 + 126], arr_index[p4*256 + 118], arr_index[p4*256 + 126], p4%2==0);
        compareSwap(arr[p4*256 + 119], arr[p4*256 + 127], arr_index[p4*256 + 119], arr_index[p4*256 + 127], p4%2==0);

        compareSwap(arr[p4*256 + 128], arr[p4*256 + 136], arr_index[p4*256 + 128], arr_index[p4*256 + 136], p4%2==0);
        compareSwap(arr[p4*256 + 129], arr[p4*256 + 137], arr_index[p4*256 + 129], arr_index[p4*256 + 137], p4%2==0);
        compareSwap(arr[p4*256 + 130], arr[p4*256 + 138], arr_index[p4*256 + 130], arr_index[p4*256 + 138], p4%2==0);
        compareSwap(arr[p4*256 + 131], arr[p4*256 + 139], arr_index[p4*256 + 131], arr_index[p4*256 + 139], p4%2==0);
        compareSwap(arr[p4*256 + 132], arr[p4*256 + 140], arr_index[p4*256 + 132], arr_index[p4*256 + 140], p4%2==0);
        compareSwap(arr[p4*256 + 133], arr[p4*256 + 141], arr_index[p4*256 + 133], arr_index[p4*256 + 141], p4%2==0);
        compareSwap(arr[p4*256 + 134], arr[p4*256 + 142], arr_index[p4*256 + 134], arr_index[p4*256 + 142], p4%2==0);
        compareSwap(arr[p4*256 + 135], arr[p4*256 + 143], arr_index[p4*256 + 135], arr_index[p4*256 + 143], p4%2==0);

        compareSwap(arr[p4*256 + 144], arr[p4*256 + 152], arr_index[p4*256 + 144], arr_index[p4*256 + 152], p4%2==0);
        compareSwap(arr[p4*256 + 145], arr[p4*256 + 153], arr_index[p4*256 + 145], arr_index[p4*256 + 153], p4%2==0);
        compareSwap(arr[p4*256 + 146], arr[p4*256 + 154], arr_index[p4*256 + 146], arr_index[p4*256 + 154], p4%2==0);
        compareSwap(arr[p4*256 + 147], arr[p4*256 + 155], arr_index[p4*256 + 147], arr_index[p4*256 + 155], p4%2==0);
        compareSwap(arr[p4*256 + 148], arr[p4*256 + 156], arr_index[p4*256 + 148], arr_index[p4*256 + 156], p4%2==0);
        compareSwap(arr[p4*256 + 149], arr[p4*256 + 157], arr_index[p4*256 + 149], arr_index[p4*256 + 157], p4%2==0);
        compareSwap(arr[p4*256 + 150], arr[p4*256 + 158], arr_index[p4*256 + 150], arr_index[p4*256 + 158], p4%2==0);
        compareSwap(arr[p4*256 + 151], arr[p4*256 + 159], arr_index[p4*256 + 151], arr_index[p4*256 + 159], p4%2==0);

        compareSwap(arr[p4*256 + 160], arr[p4*256 + 168], arr_index[p4*256 + 160], arr_index[p4*256 + 168], p4%2==0);
        compareSwap(arr[p4*256 + 161], arr[p4*256 + 169], arr_index[p4*256 + 161], arr_index[p4*256 + 169], p4%2==0);
        compareSwap(arr[p4*256 + 162], arr[p4*256 + 170], arr_index[p4*256 + 162], arr_index[p4*256 + 170], p4%2==0);
        compareSwap(arr[p4*256 + 163], arr[p4*256 + 171], arr_index[p4*256 + 163], arr_index[p4*256 + 171], p4%2==0);
        compareSwap(arr[p4*256 + 164], arr[p4*256 + 172], arr_index[p4*256 + 164], arr_index[p4*256 + 172], p4%2==0);
        compareSwap(arr[p4*256 + 165], arr[p4*256 + 173], arr_index[p4*256 + 165], arr_index[p4*256 + 173], p4%2==0);
        compareSwap(arr[p4*256 + 166], arr[p4*256 + 174], arr_index[p4*256 + 166], arr_index[p4*256 + 174], p4%2==0);
        compareSwap(arr[p4*256 + 167], arr[p4*256 + 175], arr_index[p4*256 + 167], arr_index[p4*256 + 175], p4%2==0);

        compareSwap(arr[p4*256 + 176], arr[p4*256 + 184], arr_index[p4*256 + 176], arr_index[p4*256 + 184], p4%2==0);
        compareSwap(arr[p4*256 + 177], arr[p4*256 + 185], arr_index[p4*256 + 177], arr_index[p4*256 + 185], p4%2==0);
        compareSwap(arr[p4*256 + 178], arr[p4*256 + 186], arr_index[p4*256 + 178], arr_index[p4*256 + 186], p4%2==0);
        compareSwap(arr[p4*256 + 179], arr[p4*256 + 187], arr_index[p4*256 + 179], arr_index[p4*256 + 187], p4%2==0);
        compareSwap(arr[p4*256 + 180], arr[p4*256 + 188], arr_index[p4*256 + 180], arr_index[p4*256 + 188], p4%2==0);
        compareSwap(arr[p4*256 + 181], arr[p4*256 + 189], arr_index[p4*256 + 181], arr_index[p4*256 + 189], p4%2==0);
        compareSwap(arr[p4*256 + 182], arr[p4*256 + 190], arr_index[p4*256 + 182], arr_index[p4*256 + 190], p4%2==0);
        compareSwap(arr[p4*256 + 183], arr[p4*256 + 191], arr_index[p4*256 + 183], arr_index[p4*256 + 191], p4%2==0);

        compareSwap(arr[p4*256 + 192], arr[p4*256 + 200], arr_index[p4*256 + 192], arr_index[p4*256 + 200], p4%2==0);
        compareSwap(arr[p4*256 + 193], arr[p4*256 + 201], arr_index[p4*256 + 193], arr_index[p4*256 + 201], p4%2==0);
        compareSwap(arr[p4*256 + 194], arr[p4*256 + 202], arr_index[p4*256 + 194], arr_index[p4*256 + 202], p4%2==0);
        compareSwap(arr[p4*256 + 195], arr[p4*256 + 203], arr_index[p4*256 + 195], arr_index[p4*256 + 203], p4%2==0);
        compareSwap(arr[p4*256 + 196], arr[p4*256 + 204], arr_index[p4*256 + 196], arr_index[p4*256 + 204], p4%2==0);
        compareSwap(arr[p4*256 + 197], arr[p4*256 + 205], arr_index[p4*256 + 197], arr_index[p4*256 + 205], p4%2==0);
        compareSwap(arr[p4*256 + 198], arr[p4*256 + 206], arr_index[p4*256 + 198], arr_index[p4*256 + 206], p4%2==0);
        compareSwap(arr[p4*256 + 199], arr[p4*256 + 207], arr_index[p4*256 + 199], arr_index[p4*256 + 207], p4%2==0);

        compareSwap(arr[p4*256 + 208], arr[p4*256 + 216], arr_index[p4*256 + 208], arr_index[p4*256 + 216], p4%2==0);
        compareSwap(arr[p4*256 + 209], arr[p4*256 + 217], arr_index[p4*256 + 209], arr_index[p4*256 + 217], p4%2==0);
        compareSwap(arr[p4*256 + 210], arr[p4*256 + 218], arr_index[p4*256 + 210], arr_index[p4*256 + 218], p4%2==0);
        compareSwap(arr[p4*256 + 211], arr[p4*256 + 219], arr_index[p4*256 + 211], arr_index[p4*256 + 219], p4%2==0);
        compareSwap(arr[p4*256 + 212], arr[p4*256 + 220], arr_index[p4*256 + 212], arr_index[p4*256 + 220], p4%2==0);
        compareSwap(arr[p4*256 + 213], arr[p4*256 + 221], arr_index[p4*256 + 213], arr_index[p4*256 + 221], p4%2==0);
        compareSwap(arr[p4*256 + 214], arr[p4*256 + 222], arr_index[p4*256 + 214], arr_index[p4*256 + 222], p4%2==0);
        compareSwap(arr[p4*256 + 215], arr[p4*256 + 223], arr_index[p4*256 + 215], arr_index[p4*256 + 223], p4%2==0);

        compareSwap(arr[p4*256 + 224], arr[p4*256 + 232], arr_index[p4*256 + 224], arr_index[p4*256 + 232], p4%2==0);
        compareSwap(arr[p4*256 + 225], arr[p4*256 + 233], arr_index[p4*256 + 225], arr_index[p4*256 + 233], p4%2==0);
        compareSwap(arr[p4*256 + 226], arr[p4*256 + 234], arr_index[p4*256 + 226], arr_index[p4*256 + 234], p4%2==0);
        compareSwap(arr[p4*256 + 227], arr[p4*256 + 235], arr_index[p4*256 + 227], arr_index[p4*256 + 235], p4%2==0);
        compareSwap(arr[p4*256 + 228], arr[p4*256 + 236], arr_index[p4*256 + 228], arr_index[p4*256 + 236], p4%2==0);
        compareSwap(arr[p4*256 + 229], arr[p4*256 + 237], arr_index[p4*256 + 229], arr_index[p4*256 + 237], p4%2==0);
        compareSwap(arr[p4*256 + 230], arr[p4*256 + 238], arr_index[p4*256 + 230], arr_index[p4*256 + 238], p4%2==0);
        compareSwap(arr[p4*256 + 231], arr[p4*256 + 239], arr_index[p4*256 + 231], arr_index[p4*256 + 239], p4%2==0);

        compareSwap(arr[p4*256 + 240], arr[p4*256 + 248], arr_index[p4*256 + 240], arr_index[p4*256 + 248], p4%2==0);
        compareSwap(arr[p4*256 + 241], arr[p4*256 + 249], arr_index[p4*256 + 241], arr_index[p4*256 + 249], p4%2==0);
        compareSwap(arr[p4*256 + 242], arr[p4*256 + 250], arr_index[p4*256 + 242], arr_index[p4*256 + 250], p4%2==0);
        compareSwap(arr[p4*256 + 243], arr[p4*256 + 251], arr_index[p4*256 + 243], arr_index[p4*256 + 251], p4%2==0);
        compareSwap(arr[p4*256 + 244], arr[p4*256 + 252], arr_index[p4*256 + 244], arr_index[p4*256 + 252], p4%2==0);
        compareSwap(arr[p4*256 + 245], arr[p4*256 + 253], arr_index[p4*256 + 245], arr_index[p4*256 + 253], p4%2==0);
        compareSwap(arr[p4*256 + 246], arr[p4*256 + 254], arr_index[p4*256 + 246], arr_index[p4*256 + 254], p4%2==0);
        compareSwap(arr[p4*256 + 247], arr[p4*256 + 255], arr_index[p4*256 + 247], arr_index[p4*256 + 255], p4%2==0);
    }

    for(int p5 = 0; p5 < 1; p5++)
    {
        compareSwap(arr[p5*256 + 0], arr[p5*256 + 4], arr_index[p5*256 + 0], arr_index[p5*256 + 4], p5%2==0);
        compareSwap(arr[p5*256 + 1], arr[p5*256 + 5], arr_index[p5*256 + 1], arr_index[p5*256 + 5], p5%2==0);
        compareSwap(arr[p5*256 + 2], arr[p5*256 + 6], arr_index[p5*256 + 2], arr_index[p5*256 + 6], p5%2==0);
        compareSwap(arr[p5*256 + 3], arr[p5*256 + 7], arr_index[p5*256 + 3], arr_index[p5*256 + 7], p5%2==0);

        compareSwap(arr[p5*256 + 8], arr[p5*256 + 12], arr_index[p5*256 + 8], arr_index[p5*256 + 12], p5%2==0);
        compareSwap(arr[p5*256 + 9], arr[p5*256 + 13], arr_index[p5*256 + 9], arr_index[p5*256 + 13], p5%2==0);
        compareSwap(arr[p5*256 + 10], arr[p5*256 + 14], arr_index[p5*256 + 10], arr_index[p5*256 + 14], p5%2==0);
        compareSwap(arr[p5*256 + 11], arr[p5*256 + 15], arr_index[p5*256 + 11], arr_index[p5*256 + 15], p5%2==0);

        compareSwap(arr[p5*256 + 16], arr[p5*256 + 20], arr_index[p5*256 + 16], arr_index[p5*256 + 20], p5%2==0);
        compareSwap(arr[p5*256 + 17], arr[p5*256 + 21], arr_index[p5*256 + 17], arr_index[p5*256 + 21], p5%2==0);
        compareSwap(arr[p5*256 + 18], arr[p5*256 + 22], arr_index[p5*256 + 18], arr_index[p5*256 + 22], p5%2==0);
        compareSwap(arr[p5*256 + 19], arr[p5*256 + 23], arr_index[p5*256 + 19], arr_index[p5*256 + 23], p5%2==0);

        compareSwap(arr[p5*256 + 24], arr[p5*256 + 28], arr_index[p5*256 + 24], arr_index[p5*256 + 28], p5%2==0);
        compareSwap(arr[p5*256 + 25], arr[p5*256 + 29], arr_index[p5*256 + 25], arr_index[p5*256 + 29], p5%2==0);
        compareSwap(arr[p5*256 + 26], arr[p5*256 + 30], arr_index[p5*256 + 26], arr_index[p5*256 + 30], p5%2==0);
        compareSwap(arr[p5*256 + 27], arr[p5*256 + 31], arr_index[p5*256 + 27], arr_index[p5*256 + 31], p5%2==0);

        compareSwap(arr[p5*256 + 32], arr[p5*256 + 36], arr_index[p5*256 + 32], arr_index[p5*256 + 36], p5%2==0);
        compareSwap(arr[p5*256 + 33], arr[p5*256 + 37], arr_index[p5*256 + 33], arr_index[p5*256 + 37], p5%2==0);
        compareSwap(arr[p5*256 + 34], arr[p5*256 + 38], arr_index[p5*256 + 34], arr_index[p5*256 + 38], p5%2==0);
        compareSwap(arr[p5*256 + 35], arr[p5*256 + 39], arr_index[p5*256 + 35], arr_index[p5*256 + 39], p5%2==0);

        compareSwap(arr[p5*256 + 40], arr[p5*256 + 44], arr_index[p5*256 + 40], arr_index[p5*256 + 44], p5%2==0);
        compareSwap(arr[p5*256 + 41], arr[p5*256 + 45], arr_index[p5*256 + 41], arr_index[p5*256 + 45], p5%2==0);
        compareSwap(arr[p5*256 + 42], arr[p5*256 + 46], arr_index[p5*256 + 42], arr_index[p5*256 + 46], p5%2==0);
        compareSwap(arr[p5*256 + 43], arr[p5*256 + 47], arr_index[p5*256 + 43], arr_index[p5*256 + 47], p5%2==0);

        compareSwap(arr[p5*256 + 48], arr[p5*256 + 52], arr_index[p5*256 + 48], arr_index[p5*256 + 52], p5%2==0);
        compareSwap(arr[p5*256 + 49], arr[p5*256 + 53], arr_index[p5*256 + 49], arr_index[p5*256 + 53], p5%2==0);
        compareSwap(arr[p5*256 + 50], arr[p5*256 + 54], arr_index[p5*256 + 50], arr_index[p5*256 + 54], p5%2==0);
        compareSwap(arr[p5*256 + 51], arr[p5*256 + 55], arr_index[p5*256 + 51], arr_index[p5*256 + 55], p5%2==0);

        compareSwap(arr[p5*256 + 56], arr[p5*256 + 60], arr_index[p5*256 + 56], arr_index[p5*256 + 60], p5%2==0);
        compareSwap(arr[p5*256 + 57], arr[p5*256 + 61], arr_index[p5*256 + 57], arr_index[p5*256 + 61], p5%2==0);
        compareSwap(arr[p5*256 + 58], arr[p5*256 + 62], arr_index[p5*256 + 58], arr_index[p5*256 + 62], p5%2==0);
        compareSwap(arr[p5*256 + 59], arr[p5*256 + 63], arr_index[p5*256 + 59], arr_index[p5*256 + 63], p5%2==0);

        compareSwap(arr[p5*256 + 64], arr[p5*256 + 68], arr_index[p5*256 + 64], arr_index[p5*256 + 68], p5%2==0);
        compareSwap(arr[p5*256 + 65], arr[p5*256 + 69], arr_index[p5*256 + 65], arr_index[p5*256 + 69], p5%2==0);
        compareSwap(arr[p5*256 + 66], arr[p5*256 + 70], arr_index[p5*256 + 66], arr_index[p5*256 + 70], p5%2==0);
        compareSwap(arr[p5*256 + 67], arr[p5*256 + 71], arr_index[p5*256 + 67], arr_index[p5*256 + 71], p5%2==0);

        compareSwap(arr[p5*256 + 72], arr[p5*256 + 76], arr_index[p5*256 + 72], arr_index[p5*256 + 76], p5%2==0);
        compareSwap(arr[p5*256 + 73], arr[p5*256 + 77], arr_index[p5*256 + 73], arr_index[p5*256 + 77], p5%2==0);
        compareSwap(arr[p5*256 + 74], arr[p5*256 + 78], arr_index[p5*256 + 74], arr_index[p5*256 + 78], p5%2==0);
        compareSwap(arr[p5*256 + 75], arr[p5*256 + 79], arr_index[p5*256 + 75], arr_index[p5*256 + 79], p5%2==0);

        compareSwap(arr[p5*256 + 80], arr[p5*256 + 84], arr_index[p5*256 + 80], arr_index[p5*256 + 84], p5%2==0);
        compareSwap(arr[p5*256 + 81], arr[p5*256 + 85], arr_index[p5*256 + 81], arr_index[p5*256 + 85], p5%2==0);
        compareSwap(arr[p5*256 + 82], arr[p5*256 + 86], arr_index[p5*256 + 82], arr_index[p5*256 + 86], p5%2==0);
        compareSwap(arr[p5*256 + 83], arr[p5*256 + 87], arr_index[p5*256 + 83], arr_index[p5*256 + 87], p5%2==0);

        compareSwap(arr[p5*256 + 88], arr[p5*256 + 92], arr_index[p5*256 + 88], arr_index[p5*256 + 92], p5%2==0);
        compareSwap(arr[p5*256 + 89], arr[p5*256 + 93], arr_index[p5*256 + 89], arr_index[p5*256 + 93], p5%2==0);
        compareSwap(arr[p5*256 + 90], arr[p5*256 + 94], arr_index[p5*256 + 90], arr_index[p5*256 + 94], p5%2==0);
        compareSwap(arr[p5*256 + 91], arr[p5*256 + 95], arr_index[p5*256 + 91], arr_index[p5*256 + 95], p5%2==0);

        compareSwap(arr[p5*256 + 96], arr[p5*256 + 100], arr_index[p5*256 + 96], arr_index[p5*256 + 100], p5%2==0);
        compareSwap(arr[p5*256 + 97], arr[p5*256 + 101], arr_index[p5*256 + 97], arr_index[p5*256 + 101], p5%2==0);
        compareSwap(arr[p5*256 + 98], arr[p5*256 + 102], arr_index[p5*256 + 98], arr_index[p5*256 + 102], p5%2==0);
        compareSwap(arr[p5*256 + 99], arr[p5*256 + 103], arr_index[p5*256 + 99], arr_index[p5*256 + 103], p5%2==0);

        compareSwap(arr[p5*256 + 104], arr[p5*256 + 108], arr_index[p5*256 + 104], arr_index[p5*256 + 108], p5%2==0);
        compareSwap(arr[p5*256 + 105], arr[p5*256 + 109], arr_index[p5*256 + 105], arr_index[p5*256 + 109], p5%2==0);
        compareSwap(arr[p5*256 + 106], arr[p5*256 + 110], arr_index[p5*256 + 106], arr_index[p5*256 + 110], p5%2==0);
        compareSwap(arr[p5*256 + 107], arr[p5*256 + 111], arr_index[p5*256 + 107], arr_index[p5*256 + 111], p5%2==0);

        compareSwap(arr[p5*256 + 112], arr[p5*256 + 116], arr_index[p5*256 + 112], arr_index[p5*256 + 116], p5%2==0);
        compareSwap(arr[p5*256 + 113], arr[p5*256 + 117], arr_index[p5*256 + 113], arr_index[p5*256 + 117], p5%2==0);
        compareSwap(arr[p5*256 + 114], arr[p5*256 + 118], arr_index[p5*256 + 114], arr_index[p5*256 + 118], p5%2==0);
        compareSwap(arr[p5*256 + 115], arr[p5*256 + 119], arr_index[p5*256 + 115], arr_index[p5*256 + 119], p5%2==0);

        compareSwap(arr[p5*256 + 120], arr[p5*256 + 124], arr_index[p5*256 + 120], arr_index[p5*256 + 124], p5%2==0);
        compareSwap(arr[p5*256 + 121], arr[p5*256 + 125], arr_index[p5*256 + 121], arr_index[p5*256 + 125], p5%2==0);
        compareSwap(arr[p5*256 + 122], arr[p5*256 + 126], arr_index[p5*256 + 122], arr_index[p5*256 + 126], p5%2==0);
        compareSwap(arr[p5*256 + 123], arr[p5*256 + 127], arr_index[p5*256 + 123], arr_index[p5*256 + 127], p5%2==0);

        compareSwap(arr[p5*256 + 128], arr[p5*256 + 132], arr_index[p5*256 + 128], arr_index[p5*256 + 132], p5%2==0);
        compareSwap(arr[p5*256 + 129], arr[p5*256 + 133], arr_index[p5*256 + 129], arr_index[p5*256 + 133], p5%2==0);
        compareSwap(arr[p5*256 + 130], arr[p5*256 + 134], arr_index[p5*256 + 130], arr_index[p5*256 + 134], p5%2==0);
        compareSwap(arr[p5*256 + 131], arr[p5*256 + 135], arr_index[p5*256 + 131], arr_index[p5*256 + 135], p5%2==0);

        compareSwap(arr[p5*256 + 136], arr[p5*256 + 140], arr_index[p5*256 + 136], arr_index[p5*256 + 140], p5%2==0);
        compareSwap(arr[p5*256 + 137], arr[p5*256 + 141], arr_index[p5*256 + 137], arr_index[p5*256 + 141], p5%2==0);
        compareSwap(arr[p5*256 + 138], arr[p5*256 + 142], arr_index[p5*256 + 138], arr_index[p5*256 + 142], p5%2==0);
        compareSwap(arr[p5*256 + 139], arr[p5*256 + 143], arr_index[p5*256 + 139], arr_index[p5*256 + 143], p5%2==0);

        compareSwap(arr[p5*256 + 144], arr[p5*256 + 148], arr_index[p5*256 + 144], arr_index[p5*256 + 148], p5%2==0);
        compareSwap(arr[p5*256 + 145], arr[p5*256 + 149], arr_index[p5*256 + 145], arr_index[p5*256 + 149], p5%2==0);
        compareSwap(arr[p5*256 + 146], arr[p5*256 + 150], arr_index[p5*256 + 146], arr_index[p5*256 + 150], p5%2==0);
        compareSwap(arr[p5*256 + 147], arr[p5*256 + 151], arr_index[p5*256 + 147], arr_index[p5*256 + 151], p5%2==0);

        compareSwap(arr[p5*256 + 152], arr[p5*256 + 156], arr_index[p5*256 + 152], arr_index[p5*256 + 156], p5%2==0);
        compareSwap(arr[p5*256 + 153], arr[p5*256 + 157], arr_index[p5*256 + 153], arr_index[p5*256 + 157], p5%2==0);
        compareSwap(arr[p5*256 + 154], arr[p5*256 + 158], arr_index[p5*256 + 154], arr_index[p5*256 + 158], p5%2==0);
        compareSwap(arr[p5*256 + 155], arr[p5*256 + 159], arr_index[p5*256 + 155], arr_index[p5*256 + 159], p5%2==0);

        compareSwap(arr[p5*256 + 160], arr[p5*256 + 164], arr_index[p5*256 + 160], arr_index[p5*256 + 164], p5%2==0);
        compareSwap(arr[p5*256 + 161], arr[p5*256 + 165], arr_index[p5*256 + 161], arr_index[p5*256 + 165], p5%2==0);
        compareSwap(arr[p5*256 + 162], arr[p5*256 + 166], arr_index[p5*256 + 162], arr_index[p5*256 + 166], p5%2==0);
        compareSwap(arr[p5*256 + 163], arr[p5*256 + 167], arr_index[p5*256 + 163], arr_index[p5*256 + 167], p5%2==0);

        compareSwap(arr[p5*256 + 168], arr[p5*256 + 172], arr_index[p5*256 + 168], arr_index[p5*256 + 172], p5%2==0);
        compareSwap(arr[p5*256 + 169], arr[p5*256 + 173], arr_index[p5*256 + 169], arr_index[p5*256 + 173], p5%2==0);
        compareSwap(arr[p5*256 + 170], arr[p5*256 + 174], arr_index[p5*256 + 170], arr_index[p5*256 + 174], p5%2==0);
        compareSwap(arr[p5*256 + 171], arr[p5*256 + 175], arr_index[p5*256 + 171], arr_index[p5*256 + 175], p5%2==0);

        compareSwap(arr[p5*256 + 176], arr[p5*256 + 180], arr_index[p5*256 + 176], arr_index[p5*256 + 180], p5%2==0);
        compareSwap(arr[p5*256 + 177], arr[p5*256 + 181], arr_index[p5*256 + 177], arr_index[p5*256 + 181], p5%2==0);
        compareSwap(arr[p5*256 + 178], arr[p5*256 + 182], arr_index[p5*256 + 178], arr_index[p5*256 + 182], p5%2==0);
        compareSwap(arr[p5*256 + 179], arr[p5*256 + 183], arr_index[p5*256 + 179], arr_index[p5*256 + 183], p5%2==0);

        compareSwap(arr[p5*256 + 184], arr[p5*256 + 188], arr_index[p5*256 + 184], arr_index[p5*256 + 188], p5%2==0);
        compareSwap(arr[p5*256 + 185], arr[p5*256 + 189], arr_index[p5*256 + 185], arr_index[p5*256 + 189], p5%2==0);
        compareSwap(arr[p5*256 + 186], arr[p5*256 + 190], arr_index[p5*256 + 186], arr_index[p5*256 + 190], p5%2==0);
        compareSwap(arr[p5*256 + 187], arr[p5*256 + 191], arr_index[p5*256 + 187], arr_index[p5*256 + 191], p5%2==0);

        compareSwap(arr[p5*256 + 192], arr[p5*256 + 196], arr_index[p5*256 + 192], arr_index[p5*256 + 196], p5%2==0);
        compareSwap(arr[p5*256 + 193], arr[p5*256 + 197], arr_index[p5*256 + 193], arr_index[p5*256 + 197], p5%2==0);
        compareSwap(arr[p5*256 + 194], arr[p5*256 + 198], arr_index[p5*256 + 194], arr_index[p5*256 + 198], p5%2==0);
        compareSwap(arr[p5*256 + 195], arr[p5*256 + 199], arr_index[p5*256 + 195], arr_index[p5*256 + 199], p5%2==0);

        compareSwap(arr[p5*256 + 200], arr[p5*256 + 204], arr_index[p5*256 + 200], arr_index[p5*256 + 204], p5%2==0);
        compareSwap(arr[p5*256 + 201], arr[p5*256 + 205], arr_index[p5*256 + 201], arr_index[p5*256 + 205], p5%2==0);
        compareSwap(arr[p5*256 + 202], arr[p5*256 + 206], arr_index[p5*256 + 202], arr_index[p5*256 + 206], p5%2==0);
        compareSwap(arr[p5*256 + 203], arr[p5*256 + 207], arr_index[p5*256 + 203], arr_index[p5*256 + 207], p5%2==0);

        compareSwap(arr[p5*256 + 208], arr[p5*256 + 212], arr_index[p5*256 + 208], arr_index[p5*256 + 212], p5%2==0);
        compareSwap(arr[p5*256 + 209], arr[p5*256 + 213], arr_index[p5*256 + 209], arr_index[p5*256 + 213], p5%2==0);
        compareSwap(arr[p5*256 + 210], arr[p5*256 + 214], arr_index[p5*256 + 210], arr_index[p5*256 + 214], p5%2==0);
        compareSwap(arr[p5*256 + 211], arr[p5*256 + 215], arr_index[p5*256 + 211], arr_index[p5*256 + 215], p5%2==0);

        compareSwap(arr[p5*256 + 216], arr[p5*256 + 220], arr_index[p5*256 + 216], arr_index[p5*256 + 220], p5%2==0);
        compareSwap(arr[p5*256 + 217], arr[p5*256 + 221], arr_index[p5*256 + 217], arr_index[p5*256 + 221], p5%2==0);
        compareSwap(arr[p5*256 + 218], arr[p5*256 + 222], arr_index[p5*256 + 218], arr_index[p5*256 + 222], p5%2==0);
        compareSwap(arr[p5*256 + 219], arr[p5*256 + 223], arr_index[p5*256 + 219], arr_index[p5*256 + 223], p5%2==0);

        compareSwap(arr[p5*256 + 224], arr[p5*256 + 228], arr_index[p5*256 + 224], arr_index[p5*256 + 228], p5%2==0);
        compareSwap(arr[p5*256 + 225], arr[p5*256 + 229], arr_index[p5*256 + 225], arr_index[p5*256 + 229], p5%2==0);
        compareSwap(arr[p5*256 + 226], arr[p5*256 + 230], arr_index[p5*256 + 226], arr_index[p5*256 + 230], p5%2==0);
        compareSwap(arr[p5*256 + 227], arr[p5*256 + 231], arr_index[p5*256 + 227], arr_index[p5*256 + 231], p5%2==0);

        compareSwap(arr[p5*256 + 232], arr[p5*256 + 236], arr_index[p5*256 + 232], arr_index[p5*256 + 236], p5%2==0);
        compareSwap(arr[p5*256 + 233], arr[p5*256 + 237], arr_index[p5*256 + 233], arr_index[p5*256 + 237], p5%2==0);
        compareSwap(arr[p5*256 + 234], arr[p5*256 + 238], arr_index[p5*256 + 234], arr_index[p5*256 + 238], p5%2==0);
        compareSwap(arr[p5*256 + 235], arr[p5*256 + 239], arr_index[p5*256 + 235], arr_index[p5*256 + 239], p5%2==0);

        compareSwap(arr[p5*256 + 240], arr[p5*256 + 244], arr_index[p5*256 + 240], arr_index[p5*256 + 244], p5%2==0);
        compareSwap(arr[p5*256 + 241], arr[p5*256 + 245], arr_index[p5*256 + 241], arr_index[p5*256 + 245], p5%2==0);
        compareSwap(arr[p5*256 + 242], arr[p5*256 + 246], arr_index[p5*256 + 242], arr_index[p5*256 + 246], p5%2==0);
        compareSwap(arr[p5*256 + 243], arr[p5*256 + 247], arr_index[p5*256 + 243], arr_index[p5*256 + 247], p5%2==0);

        compareSwap(arr[p5*256 + 248], arr[p5*256 + 252], arr_index[p5*256 + 248], arr_index[p5*256 + 252], p5%2==0);
        compareSwap(arr[p5*256 + 249], arr[p5*256 + 253], arr_index[p5*256 + 249], arr_index[p5*256 + 253], p5%2==0);
        compareSwap(arr[p5*256 + 250], arr[p5*256 + 254], arr_index[p5*256 + 250], arr_index[p5*256 + 254], p5%2==0);
        compareSwap(arr[p5*256 + 251], arr[p5*256 + 255], arr_index[p5*256 + 251], arr_index[p5*256 + 255], p5%2==0);
    }

    for(int p6 = 0; p6 < 1; p6++)
    {
        compareSwap(arr[p6*256 + 0], arr[p6*256 + 2], arr_index[p6*256 + 0], arr_index[p6*256 + 2], p6%2==0);
        compareSwap(arr[p6*256 + 1], arr[p6*256 + 3], arr_index[p6*256 + 1], arr_index[p6*256 + 3], p6%2==0);

        compareSwap(arr[p6*256 + 4], arr[p6*256 + 6], arr_index[p6*256 + 4], arr_index[p6*256 + 6], p6%2==0);
        compareSwap(arr[p6*256 + 5], arr[p6*256 + 7], arr_index[p6*256 + 5], arr_index[p6*256 + 7], p6%2==0);

        compareSwap(arr[p6*256 + 8], arr[p6*256 + 10], arr_index[p6*256 + 8], arr_index[p6*256 + 10], p6%2==0);
        compareSwap(arr[p6*256 + 9], arr[p6*256 + 11], arr_index[p6*256 + 9], arr_index[p6*256 + 11], p6%2==0);

        compareSwap(arr[p6*256 + 12], arr[p6*256 + 14], arr_index[p6*256 + 12], arr_index[p6*256 + 14], p6%2==0);
        compareSwap(arr[p6*256 + 13], arr[p6*256 + 15], arr_index[p6*256 + 13], arr_index[p6*256 + 15], p6%2==0);

        compareSwap(arr[p6*256 + 16], arr[p6*256 + 18], arr_index[p6*256 + 16], arr_index[p6*256 + 18], p6%2==0);
        compareSwap(arr[p6*256 + 17], arr[p6*256 + 19], arr_index[p6*256 + 17], arr_index[p6*256 + 19], p6%2==0);

        compareSwap(arr[p6*256 + 20], arr[p6*256 + 22], arr_index[p6*256 + 20], arr_index[p6*256 + 22], p6%2==0);
        compareSwap(arr[p6*256 + 21], arr[p6*256 + 23], arr_index[p6*256 + 21], arr_index[p6*256 + 23], p6%2==0);

        compareSwap(arr[p6*256 + 24], arr[p6*256 + 26], arr_index[p6*256 + 24], arr_index[p6*256 + 26], p6%2==0);
        compareSwap(arr[p6*256 + 25], arr[p6*256 + 27], arr_index[p6*256 + 25], arr_index[p6*256 + 27], p6%2==0);
        compareSwap(arr[p6*256 + 28], arr[p6*256 + 30], arr_index[p6*256 + 28], arr_index[p6*256 + 30], p6%2==0);
        compareSwap(arr[p6*256 + 29], arr[p6*256 + 31], arr_index[p6*256 + 29], arr_index[p6*256 + 31], p6%2==0);

        compareSwap(arr[p6*256 + 32], arr[p6*256 + 34], arr_index[p6*256 + 32], arr_index[p6*256 + 34], p6%2==0);
        compareSwap(arr[p6*256 + 33], arr[p6*256 + 35], arr_index[p6*256 + 33], arr_index[p6*256 + 35], p6%2==0);

        compareSwap(arr[p6*256 + 36], arr[p6*256 + 38], arr_index[p6*256 + 36], arr_index[p6*256 + 38], p6%2==0);
        compareSwap(arr[p6*256 + 37], arr[p6*256 + 39], arr_index[p6*256 + 37], arr_index[p6*256 + 39], p6%2==0);

        compareSwap(arr[p6*256 + 40], arr[p6*256 + 42], arr_index[p6*256 + 40], arr_index[p6*256 + 42], p6%2==0);
        compareSwap(arr[p6*256 + 41], arr[p6*256 + 43], arr_index[p6*256 + 41], arr_index[p6*256 + 43], p6%2==0);
        compareSwap(arr[p6*256 + 44], arr[p6*256 + 46], arr_index[p6*256 + 44], arr_index[p6*256 + 46], p6%2==0);
        compareSwap(arr[p6*256 + 45], arr[p6*256 + 47], arr_index[p6*256 + 45], arr_index[p6*256 + 47], p6%2==0);

        compareSwap(arr[p6*256 + 48], arr[p6*256 + 50], arr_index[p6*256 + 48], arr_index[p6*256 + 50], p6%2==0);
        compareSwap(arr[p6*256 + 49], arr[p6*256 + 51], arr_index[p6*256 + 49], arr_index[p6*256 + 51], p6%2==0);

        compareSwap(arr[p6*256 + 52], arr[p6*256 + 54], arr_index[p6*256 + 52], arr_index[p6*256 + 54], p6%2==0);
        compareSwap(arr[p6*256 + 53], arr[p6*256 + 55], arr_index[p6*256 + 53], arr_index[p6*256 + 55], p6%2==0);
        compareSwap(arr[p6*256 + 56], arr[p6*256 + 58], arr_index[p6*256 + 56], arr_index[p6*256 + 58], p6%2==0);
        compareSwap(arr[p6*256 + 57], arr[p6*256 + 59], arr_index[p6*256 + 57], arr_index[p6*256 + 59], p6%2==0);

        compareSwap(arr[p6*256 + 60], arr[p6*256 + 62], arr_index[p6*256 + 60], arr_index[p6*256 + 62], p6%2==0);
        compareSwap(arr[p6*256 + 61], arr[p6*256 + 63], arr_index[p6*256 + 61], arr_index[p6*256 + 63], p6%2==0);

        compareSwap(arr[p6*256 + 64], arr[p6*256 + 66], arr_index[p6*256 + 64], arr_index[p6*256 + 66], p6%2==0);
        compareSwap(arr[p6*256 + 65], arr[p6*256 + 67], arr_index[p6*256 + 65], arr_index[p6*256 + 67], p6%2==0);

        compareSwap(arr[p6*256 + 68], arr[p6*256 + 70], arr_index[p6*256 + 68], arr_index[p6*256 + 70], p6%2==0);
        compareSwap(arr[p6*256 + 69], arr[p6*256 + 71], arr_index[p6*256 + 69], arr_index[p6*256 + 71], p6%2==0);
        compareSwap(arr[p6*256 + 72], arr[p6*256 + 74], arr_index[p6*256 + 72], arr_index[p6*256 + 74], p6%2==0);
        compareSwap(arr[p6*256 + 73], arr[p6*256 + 75], arr_index[p6*256 + 73], arr_index[p6*256 + 75], p6%2==0);

        compareSwap(arr[p6*256 + 76], arr[p6*256 + 78], arr_index[p6*256 + 76], arr_index[p6*256 + 78], p6%2==0);
        compareSwap(arr[p6*256 + 77], arr[p6*256 + 79], arr_index[p6*256 + 77], arr_index[p6*256 + 79], p6%2==0);

        compareSwap(arr[p6*256 + 80], arr[p6*256 + 82], arr_index[p6*256 + 80], arr_index[p6*256 + 82], p6%2==0);
        compareSwap(arr[p6*256 + 81], arr[p6*256 + 83], arr_index[p6*256 + 81], arr_index[p6*256 + 83], p6%2==0);

        compareSwap(arr[p6*256 + 84], arr[p6*256 + 86], arr_index[p6*256 + 84], arr_index[p6*256 + 86], p6%2==0);
        compareSwap(arr[p6*256 + 85], arr[p6*256 + 87], arr_index[p6*256 + 85], arr_index[p6*256 + 87], p6%2==0);
        compareSwap(arr[p6*256 + 88], arr[p6*256 + 90], arr_index[p6*256 + 88], arr_index[p6*256 + 90], p6%2==0);
        compareSwap(arr[p6*256 + 89], arr[p6*256 + 91], arr_index[p6*256 + 89], arr_index[p6*256 + 91], p6%2==0);

        compareSwap(arr[p6*256 + 92], arr[p6*256 + 94], arr_index[p6*256 + 92], arr_index[p6*256 + 94], p6%2==0);
        compareSwap(arr[p6*256 + 93], arr[p6*256 + 95], arr_index[p6*256 + 93], arr_index[p6*256 + 95], p6%2==0);

        compareSwap(arr[p6*256 + 96], arr[p6*256 + 98], arr_index[p6*256 + 96], arr_index[p6*256 + 98], p6%2==0);
        compareSwap(arr[p6*256 + 97], arr[p6*256 + 99], arr_index[p6*256 + 97], arr_index[p6*256 + 99], p6%2==0);

        compareSwap(arr[p6*256 + 100], arr[p6*256 + 102], arr_index[p6*256 + 100], arr_index[p6*256 + 102], p6%2==0);
        compareSwap(arr[p6*256 + 101], arr[p6*256 + 103], arr_index[p6*256 + 101], arr_index[p6*256 + 103], p6%2==0);
        compareSwap(arr[p6*256 + 104], arr[p6*256 + 106], arr_index[p6*256 + 104], arr_index[p6*256 + 106], p6%2==0);
        compareSwap(arr[p6*256 + 105], arr[p6*256 + 107], arr_index[p6*256 + 105], arr_index[p6*256 + 107], p6%2==0);

        compareSwap(arr[p6*256 + 108], arr[p6*256 + 110], arr_index[p6*256 + 108], arr_index[p6*256 + 110], p6%2==0);
        compareSwap(arr[p6*256 + 109], arr[p6*256 + 111], arr_index[p6*256 + 109], arr_index[p6*256 + 111], p6%2==0);

        compareSwap(arr[p6*256 + 112], arr[p6*256 + 114], arr_index[p6*256 + 112], arr_index[p6*256 + 114], p6%2==0);
        compareSwap(arr[p6*256 + 113], arr[p6*256 + 115], arr_index[p6*256 + 113], arr_index[p6*256 + 115], p6%2==0);

        compareSwap(arr[p6*256 + 116], arr[p6*256 + 118], arr_index[p6*256 + 116], arr_index[p6*256 + 118], p6%2==0);
        compareSwap(arr[p6*256 + 117], arr[p6*256 + 119], arr_index[p6*256 + 117], arr_index[p6*256 + 119], p6%2==0);
        compareSwap(arr[p6*256 + 120], arr[p6*256 + 122], arr_index[p6*256 + 120], arr_index[p6*256 + 122], p6%2==0);
        compareSwap(arr[p6*256 + 121], arr[p6*256 + 123], arr_index[p6*256 + 121], arr_index[p6*256 + 123], p6%2==0);

        compareSwap(arr[p6*256 + 124], arr[p6*256 + 126], arr_index[p6*256 + 124], arr_index[p6*256 + 126], p6%2==0);
        compareSwap(arr[p6*256 + 125], arr[p6*256 + 127], arr_index[p6*256 + 125], arr_index[p6*256 + 127], p6%2==0);

        compareSwap(arr[p6*256 + 128], arr[p6*256 + 130], arr_index[p6*256 + 128], arr_index[p6*256 + 130], p6%2==0);
        compareSwap(arr[p6*256 + 129], arr[p6*256 + 131], arr_index[p6*256 + 129], arr_index[p6*256 + 131], p6%2==0);

        compareSwap(arr[p6*256 + 132], arr[p6*256 + 134], arr_index[p6*256 + 132], arr_index[p6*256 + 134], p6%2==0);
        compareSwap(arr[p6*256 + 133], arr[p6*256 + 135], arr_index[p6*256 + 133], arr_index[p6*256 + 135], p6%2==0);
        compareSwap(arr[p6*256 + 136], arr[p6*256 + 138], arr_index[p6*256 + 136], arr_index[p6*256 + 138], p6%2==0);
        compareSwap(arr[p6*256 + 137], arr[p6*256 + 139], arr_index[p6*256 + 137], arr_index[p6*256 + 139], p6%2==0);

        compareSwap(arr[p6*256 + 140], arr[p6*256 + 142], arr_index[p6*256 + 140], arr_index[p6*256 + 142], p6%2==0);
        compareSwap(arr[p6*256 + 141], arr[p6*256 + 143], arr_index[p6*256 + 141], arr_index[p6*256 + 143], p6%2==0);
        compareSwap(arr[p6*256 + 144], arr[p6*256 + 146], arr_index[p6*256 + 144], arr_index[p6*256 + 146], p6%2==0);
        compareSwap(arr[p6*256 + 145], arr[p6*256 + 147], arr_index[p6*256 + 145], arr_index[p6*256 + 147], p6%2==0);

        compareSwap(arr[p6*256 + 148], arr[p6*256 + 150], arr_index[p6*256 + 148], arr_index[p6*256 + 150], p6%2==0);
        compareSwap(arr[p6*256 + 149], arr[p6*256 + 151], arr_index[p6*256 + 149], arr_index[p6*256 + 151], p6%2==0);
        compareSwap(arr[p6*256 + 152], arr[p6*256 + 154], arr_index[p6*256 + 152], arr_index[p6*256 + 154], p6%2==0);
        compareSwap(arr[p6*256 + 153], arr[p6*256 + 155], arr_index[p6*256 + 153], arr_index[p6*256 + 155], p6%2==0);

        compareSwap(arr[p6*256 + 156], arr[p6*256 + 158], arr_index[p6*256 + 156], arr_index[p6*256 + 158], p6%2==0);
        compareSwap(arr[p6*256 + 157], arr[p6*256 + 159], arr_index[p6*256 + 157], arr_index[p6*256 + 159], p6%2==0);
        compareSwap(arr[p6*256 + 160], arr[p6*256 + 162], arr_index[p6*256 + 160], arr_index[p6*256 + 162], p6%2==0);
        compareSwap(arr[p6*256 + 161], arr[p6*256 + 163], arr_index[p6*256 + 161], arr_index[p6*256 + 163], p6%2==0);

        compareSwap(arr[p6*256 + 164], arr[p6*256 + 166], arr_index[p6*256 + 164], arr_index[p6*256 + 166], p6%2==0);
        compareSwap(arr[p6*256 + 165], arr[p6*256 + 167], arr_index[p6*256 + 165], arr_index[p6*256 + 167], p6%2==0);
        compareSwap(arr[p6*256 + 168], arr[p6*256 + 170], arr_index[p6*256 + 168], arr_index[p6*256 + 170], p6%2==0);
        compareSwap(arr[p6*256 + 169], arr[p6*256 + 171], arr_index[p6*256 + 169], arr_index[p6*256 + 171], p6%2==0);

        compareSwap(arr[p6*256 + 172], arr[p6*256 + 174], arr_index[p6*256 + 172], arr_index[p6*256 + 174], p6%2==0);
        compareSwap(arr[p6*256 + 173], arr[p6*256 + 175], arr_index[p6*256 + 173], arr_index[p6*256 + 175], p6%2==0);
        compareSwap(arr[p6*256 + 176], arr[p6*256 + 178], arr_index[p6*256 + 176], arr_index[p6*256 + 178], p6%2==0);
        compareSwap(arr[p6*256 + 177], arr[p6*256 + 179], arr_index[p6*256 + 177], arr_index[p6*256 + 179], p6%2==0);

        compareSwap(arr[p6*256 + 180], arr[p6*256 + 182], arr_index[p6*256 + 180], arr_index[p6*256 + 182], p6%2==0);
        compareSwap(arr[p6*256 + 181], arr[p6*256 + 183], arr_index[p6*256 + 181], arr_index[p6*256 + 183], p6%2==0);
        compareSwap(arr[p6*256 + 184], arr[p6*256 + 186], arr_index[p6*256 + 184], arr_index[p6*256 + 186], p6%2==0);
        compareSwap(arr[p6*256 + 185], arr[p6*256 + 187], arr_index[p6*256 + 185], arr_index[p6*256 + 187], p6%2==0);

        compareSwap(arr[p6*256 + 188], arr[p6*256 + 190], arr_index[p6*256 + 188], arr_index[p6*256 + 190], p6%2==0);
        compareSwap(arr[p6*256 + 189], arr[p6*256 + 191], arr_index[p6*256 + 189], arr_index[p6*256 + 191], p6%2==0);

        compareSwap(arr[p6*256 + 192], arr[p6*256 + 194], arr_index[p6*256 + 192], arr_index[p6*256 + 194], p6%2==0);
        compareSwap(arr[p6*256 + 193], arr[p6*256 + 195], arr_index[p6*256 + 193], arr_index[p6*256 + 195], p6%2==0);

        compareSwap(arr[p6*256 + 196], arr[p6*256 + 198], arr_index[p6*256 + 196], arr_index[p6*256 + 198], p6%2==0);
        compareSwap(arr[p6*256 + 197], arr[p6*256 + 199], arr_index[p6*256 + 197], arr_index[p6*256 + 199], p6%2==0);
        compareSwap(arr[p6*256 + 200], arr[p6*256 + 202], arr_index[p6*256 + 200], arr_index[p6*256 + 202], p6%2==0);
        compareSwap(arr[p6*256 + 201], arr[p6*256 + 203], arr_index[p6*256 + 201], arr_index[p6*256 + 203], p6%2==0);

        compareSwap(arr[p6*256 + 204], arr[p6*256 + 206], arr_index[p6*256 + 204], arr_index[p6*256 + 206], p6%2==0);
        compareSwap(arr[p6*256 + 205], arr[p6*256 + 207], arr_index[p6*256 + 205], arr_index[p6*256 + 207], p6%2==0);

        compareSwap(arr[p6*256 + 208], arr[p6*256 + 210], arr_index[p6*256 + 208], arr_index[p6*256 + 210], p6%2==0);
        compareSwap(arr[p6*256 + 209], arr[p6*256 + 211], arr_index[p6*256 + 209], arr_index[p6*256 + 211], p6%2==0);

        compareSwap(arr[p6*256 + 212], arr[p6*256 + 214], arr_index[p6*256 + 212], arr_index[p6*256 + 214], p6%2==0);
        compareSwap(arr[p6*256 + 213], arr[p6*256 + 215], arr_index[p6*256 + 213], arr_index[p6*256 + 215], p6%2==0);

        compareSwap(arr[p6*256 + 216], arr[p6*256 + 218], arr_index[p6*256 + 216], arr_index[p6*256 + 218], p6%2==0);
        compareSwap(arr[p6*256 + 217], arr[p6*256 + 219], arr_index[p6*256 + 217], arr_index[p6*256 + 219], p6%2==0);
        compareSwap(arr[p6*256 + 220], arr[p6*256 + 222], arr_index[p6*256 + 220], arr_index[p6*256 + 222], p6%2==0);
        compareSwap(arr[p6*256 + 221], arr[p6*256 + 223], arr_index[p6*256 + 221], arr_index[p6*256 + 223], p6%2==0);

        compareSwap(arr[p6*256 + 224], arr[p6*256 + 226], arr_index[p6*256 + 224], arr_index[p6*256 + 226], p6%2==0);
        compareSwap(arr[p6*256 + 225], arr[p6*256 + 227], arr_index[p6*256 + 225], arr_index[p6*256 + 227], p6%2==0);
        compareSwap(arr[p6*256 + 228], arr[p6*256 + 230], arr_index[p6*256 + 228], arr_index[p6*256 + 230], p6%2==0);
        compareSwap(arr[p6*256 + 229], arr[p6*256 + 231], arr_index[p6*256 + 229], arr_index[p6*256 + 231], p6%2==0);

        compareSwap(arr[p6*256 + 232], arr[p6*256 + 234], arr_index[p6*256 + 232], arr_index[p6*256 + 234], p6%2==0);
        compareSwap(arr[p6*256 + 233], arr[p6*256 + 235], arr_index[p6*256 + 233], arr_index[p6*256 + 235], p6%2==0);
        compareSwap(arr[p6*256 + 236], arr[p6*256 + 238], arr_index[p6*256 + 236], arr_index[p6*256 + 238], p6%2==0);
        compareSwap(arr[p6*256 + 237], arr[p6*256 + 239], arr_index[p6*256 + 237], arr_index[p6*256 + 239], p6%2==0);

        compareSwap(arr[p6*256 + 240], arr[p6*256 + 242], arr_index[p6*256 + 240], arr_index[p6*256 + 242], p6%2==0);
        compareSwap(arr[p6*256 + 241], arr[p6*256 + 243], arr_index[p6*256 + 241], arr_index[p6*256 + 243], p6%2==0);
        compareSwap(arr[p6*256 + 244], arr[p6*256 + 246], arr_index[p6*256 + 244], arr_index[p6*256 + 246], p6%2==0);
        compareSwap(arr[p6*256 + 245], arr[p6*256 + 247], arr_index[p6*256 + 245], arr_index[p6*256 + 247], p6%2==0);

        compareSwap(arr[p6*256 + 248], arr[p6*256 + 250], arr_index[p6*256 + 248], arr_index[p6*256 + 250], p6%2==0);
        compareSwap(arr[p6*256 + 249], arr[p6*256 + 251], arr_index[p6*256 + 249], arr_index[p6*256 + 251], p6%2==0);
        compareSwap(arr[p6*256 + 252], arr[p6*256 + 254], arr_index[p6*256 + 252], arr_index[p6*256 + 254], p6%2==0);
        compareSwap(arr[p6*256 + 253], arr[p6*256 + 255], arr_index[p6*256 + 253], arr_index[p6*256 + 255], p6%2==0);
    }

    for(int p7 = 0; p7 < 1; p7++){
        compareSwap(arr[p7*256 + 0], arr[p7*256 + 1], arr_index[p7*256 + 0], arr_index[p7*256 + 1], p7%2==0);
        compareSwap(arr[p7*256 + 2], arr[p7*256 + 3], arr_index[p7*256 + 2], arr_index[p7*256 + 3], p7%2==0);
        compareSwap(arr[p7*256 + 4], arr[p7*256 + 5], arr_index[p7*256 + 4], arr_index[p7*256 + 5], p7%2==0);
        compareSwap(arr[p7*256 + 6], arr[p7*256 + 7], arr_index[p7*256 + 6], arr_index[p7*256 + 7], p7%2==0);
        compareSwap(arr[p7*256 + 8], arr[p7*256 + 9], arr_index[p7*256 + 8], arr_index[p7*256 + 9], p7%2==0);
        compareSwap(arr[p7*256 + 10], arr[p7*256 + 11], arr_index[p7*256 + 10], arr_index[p7*256 + 11], p7%2==0);
        compareSwap(arr[p7*256 + 12], arr[p7*256 + 13], arr_index[p7*256 + 12], arr_index[p7*256 + 13], p7%2==0);
        compareSwap(arr[p7*256 + 14], arr[p7*256 + 15], arr_index[p7*256 + 14], arr_index[p7*256 + 15], p7%2==0);
        compareSwap(arr[p7*256 + 16], arr[p7*256 + 17], arr_index[p7*256 + 16], arr_index[p7*256 + 17], p7%2==0);
        compareSwap(arr[p7*256 + 18], arr[p7*256 + 19], arr_index[p7*256 + 18], arr_index[p7*256 + 19], p7%2==0);
        compareSwap(arr[p7*256 + 20], arr[p7*256 + 21], arr_index[p7*256 + 20], arr_index[p7*256 + 21], p7%2==0);
        compareSwap(arr[p7*256 + 22], arr[p7*256 + 23], arr_index[p7*256 + 22], arr_index[p7*256 + 23], p7%2==0);
        compareSwap(arr[p7*256 + 24], arr[p7*256 + 25], arr_index[p7*256 + 24], arr_index[p7*256 + 25], p7%2==0);
        compareSwap(arr[p7*256 + 26], arr[p7*256 + 27], arr_index[p7*256 + 26], arr_index[p7*256 + 27], p7%2==0);
        compareSwap(arr[p7*256 + 28], arr[p7*256 + 29], arr_index[p7*256 + 28], arr_index[p7*256 + 29], p7%2==0);
        compareSwap(arr[p7*256 + 30], arr[p7*256 + 31], arr_index[p7*256 + 30], arr_index[p7*256 + 31], p7%2==0);
        compareSwap(arr[p7*256 + 32], arr[p7*256 + 33], arr_index[p7*256 + 32], arr_index[p7*256 + 33], p7%2==0);
        compareSwap(arr[p7*256 + 34], arr[p7*256 + 35], arr_index[p7*256 + 34], arr_index[p7*256 + 35], p7%2==0);
        compareSwap(arr[p7*256 + 36], arr[p7*256 + 37], arr_index[p7*256 + 36], arr_index[p7*256 + 37], p7%2==0);
        compareSwap(arr[p7*256 + 38], arr[p7*256 + 39], arr_index[p7*256 + 38], arr_index[p7*256 + 39], p7%2==0);
        compareSwap(arr[p7*256 + 40], arr[p7*256 + 41], arr_index[p7*256 + 40], arr_index[p7*256 + 41], p7%2==0);
        compareSwap(arr[p7*256 + 42], arr[p7*256 + 43], arr_index[p7*256 + 42], arr_index[p7*256 + 43], p7%2==0);
        compareSwap(arr[p7*256 + 44], arr[p7*256 + 45], arr_index[p7*256 + 44], arr_index[p7*256 + 45], p7%2==0);
        compareSwap(arr[p7*256 + 46], arr[p7*256 + 47], arr_index[p7*256 + 46], arr_index[p7*256 + 47], p7%2==0);
        compareSwap(arr[p7*256 + 48], arr[p7*256 + 49], arr_index[p7*256 + 48], arr_index[p7*256 + 49], p7%2==0);
        compareSwap(arr[p7*256 + 50], arr[p7*256 + 51], arr_index[p7*256 + 50], arr_index[p7*256 + 51], p7%2==0);
        compareSwap(arr[p7*256 + 52], arr[p7*256 + 53], arr_index[p7*256 + 52], arr_index[p7*256 + 53], p7%2==0);
        compareSwap(arr[p7*256 + 54], arr[p7*256 + 55], arr_index[p7*256 + 54], arr_index[p7*256 + 55], p7%2==0);
        compareSwap(arr[p7*256 + 56], arr[p7*256 + 57], arr_index[p7*256 + 56], arr_index[p7*256 + 57], p7%2==0);
        compareSwap(arr[p7*256 + 58], arr[p7*256 + 59], arr_index[p7*256 + 58], arr_index[p7*256 + 59], p7%2==0);
        compareSwap(arr[p7*256 + 60], arr[p7*256 + 61], arr_index[p7*256 + 60], arr_index[p7*256 + 61], p7%2==0);
        compareSwap(arr[p7*256 + 62], arr[p7*256 + 63], arr_index[p7*256 + 62], arr_index[p7*256 + 63], p7%2==0);

        compareSwap(arr[p7*256 + 64], arr[p7*256 + 65], arr_index[p7*256 + 64], arr_index[p7*256 + 65], p7%2==0);
        compareSwap(arr[p7*256 + 66], arr[p7*256 + 67], arr_index[p7*256 + 66], arr_index[p7*256 + 67], p7%2==0);
        compareSwap(arr[p7*256 + 68], arr[p7*256 + 69], arr_index[p7*256 + 68], arr_index[p7*256 + 69], p7%2==0);
        compareSwap(arr[p7*256 + 70], arr[p7*256 + 71], arr_index[p7*256 + 70], arr_index[p7*256 + 71], p7%2==0);
        compareSwap(arr[p7*256 + 72], arr[p7*256 + 73], arr_index[p7*256 + 72], arr_index[p7*256 + 73], p7%2==0);
        compareSwap(arr[p7*256 + 74], arr[p7*256 + 75], arr_index[p7*256 + 74], arr_index[p7*256 + 75], p7%2==0);
        compareSwap(arr[p7*256 + 76], arr[p7*256 + 77], arr_index[p7*256 + 76], arr_index[p7*256 + 77], p7%2==0);
        compareSwap(arr[p7*256 + 78], arr[p7*256 + 79], arr_index[p7*256 + 78], arr_index[p7*256 + 79], p7%2==0);
        compareSwap(arr[p7*256 + 80], arr[p7*256 + 81], arr_index[p7*256 + 80], arr_index[p7*256 + 81], p7%2==0);
        compareSwap(arr[p7*256 + 82], arr[p7*256 + 83], arr_index[p7*256 + 82], arr_index[p7*256 + 83], p7%2==0);
        compareSwap(arr[p7*256 + 84], arr[p7*256 + 85], arr_index[p7*256 + 84], arr_index[p7*256 + 85], p7%2==0);
        compareSwap(arr[p7*256 + 86], arr[p7*256 + 87], arr_index[p7*256 + 86], arr_index[p7*256 + 87], p7%2==0);
        compareSwap(arr[p7*256 + 88], arr[p7*256 + 89], arr_index[p7*256 + 88], arr_index[p7*256 + 89], p7%2==0);
        compareSwap(arr[p7*256 + 90], arr[p7*256 + 91], arr_index[p7*256 + 90], arr_index[p7*256 + 91], p7%2==0);
        compareSwap(arr[p7*256 + 92], arr[p7*256 + 93], arr_index[p7*256 + 92], arr_index[p7*256 + 93], p7%2==0);
        compareSwap(arr[p7*256 + 94], arr[p7*256 + 95], arr_index[p7*256 + 94], arr_index[p7*256 + 95], p7%2==0);
        compareSwap(arr[p7*256 + 96], arr[p7*256 + 97], arr_index[p7*256 + 96], arr_index[p7*256 + 97], p7%2==0);
        compareSwap(arr[p7*256 + 98], arr[p7*256 + 99], arr_index[p7*256 + 98], arr_index[p7*256 + 99], p7%2==0);
        compareSwap(arr[p7*256 + 100], arr[p7*256 + 101], arr_index[p7*256 + 100], arr_index[p7*256 + 101], p7%2==0);
        compareSwap(arr[p7*256 + 102], arr[p7*256 + 103], arr_index[p7*256 + 102], arr_index[p7*256 + 103], p7%2==0);
        compareSwap(arr[p7*256 + 104], arr[p7*256 + 105], arr_index[p7*256 + 104], arr_index[p7*256 + 105], p7%2==0);
        compareSwap(arr[p7*256 + 106], arr[p7*256 + 107], arr_index[p7*256 + 106], arr_index[p7*256 + 107], p7%2==0);
        compareSwap(arr[p7*256 + 108], arr[p7*256 + 109], arr_index[p7*256 + 108], arr_index[p7*256 + 109], p7%2==0);
        compareSwap(arr[p7*256 + 110], arr[p7*256 + 111], arr_index[p7*256 + 110], arr_index[p7*256 + 111], p7%2==0);
        compareSwap(arr[p7*256 + 112], arr[p7*256 + 113], arr_index[p7*256 + 112], arr_index[p7*256 + 113], p7%2==0);
        compareSwap(arr[p7*256 + 114], arr[p7*256 + 115], arr_index[p7*256 + 114], arr_index[p7*256 + 115], p7%2==0);
        compareSwap(arr[p7*256 + 116], arr[p7*256 + 117], arr_index[p7*256 + 116], arr_index[p7*256 + 117], p7%2==0);
        compareSwap(arr[p7*256 + 118], arr[p7*256 + 119], arr_index[p7*256 + 118], arr_index[p7*256 + 119], p7%2==0);
        compareSwap(arr[p7*256 + 120], arr[p7*256 + 121], arr_index[p7*256 + 120], arr_index[p7*256 + 121], p7%2==0);
        compareSwap(arr[p7*256 + 122], arr[p7*256 + 123], arr_index[p7*256 + 122], arr_index[p7*256 + 123], p7%2==0);
        compareSwap(arr[p7*256 + 124], arr[p7*256 + 125], arr_index[p7*256 + 124], arr_index[p7*256 + 125], p7%2==0);
        compareSwap(arr[p7*256 + 126], arr[p7*256 + 127], arr_index[p7*256 + 126], arr_index[p7*256 + 127], p7%2==0);

        compareSwap(arr[p7*256 + 128], arr[p7*256 + 129], arr_index[p7*256 + 128], arr_index[p7*256 + 129], p7%2==0);
        compareSwap(arr[p7*256 + 130], arr[p7*256 + 131], arr_index[p7*256 + 130], arr_index[p7*256 + 131], p7%2==0);
        compareSwap(arr[p7*256 + 132], arr[p7*256 + 133], arr_index[p7*256 + 132], arr_index[p7*256 + 133], p7%2==0);
        compareSwap(arr[p7*256 + 134], arr[p7*256 + 135], arr_index[p7*256 + 134], arr_index[p7*256 + 135], p7%2==0);
        compareSwap(arr[p7*256 + 136], arr[p7*256 + 137], arr_index[p7*256 + 136], arr_index[p7*256 + 137], p7%2==0);
        compareSwap(arr[p7*256 + 138], arr[p7*256 + 139], arr_index[p7*256 + 138], arr_index[p7*256 + 139], p7%2==0);
        compareSwap(arr[p7*256 + 140], arr[p7*256 + 141], arr_index[p7*256 + 140], arr_index[p7*256 + 141], p7%2==0);
        compareSwap(arr[p7*256 + 142], arr[p7*256 + 143], arr_index[p7*256 + 142], arr_index[p7*256 + 143], p7%2==0);
        compareSwap(arr[p7*256 + 144], arr[p7*256 + 145], arr_index[p7*256 + 144], arr_index[p7*256 + 145], p7%2==0);
        compareSwap(arr[p7*256 + 146], arr[p7*256 + 147], arr_index[p7*256 + 146], arr_index[p7*256 + 147], p7%2==0);
        compareSwap(arr[p7*256 + 148], arr[p7*256 + 149], arr_index[p7*256 + 148], arr_index[p7*256 + 149], p7%2==0);
        compareSwap(arr[p7*256 + 150], arr[p7*256 + 151], arr_index[p7*256 + 150], arr_index[p7*256 + 151], p7%2==0);
        compareSwap(arr[p7*256 + 152], arr[p7*256 + 153], arr_index[p7*256 + 152], arr_index[p7*256 + 153], p7%2==0);
        compareSwap(arr[p7*256 + 154], arr[p7*256 + 155], arr_index[p7*256 + 154], arr_index[p7*256 + 155], p7%2==0);
        compareSwap(arr[p7*256 + 156], arr[p7*256 + 157], arr_index[p7*256 + 156], arr_index[p7*256 + 157], p7%2==0);
        compareSwap(arr[p7*256 + 158], arr[p7*256 + 159], arr_index[p7*256 + 158], arr_index[p7*256 + 159], p7%2==0);
        compareSwap(arr[p7*256 + 160], arr[p7*256 + 161], arr_index[p7*256 + 160], arr_index[p7*256 + 161], p7%2==0);
        compareSwap(arr[p7*256 + 162], arr[p7*256 + 163], arr_index[p7*256 + 162], arr_index[p7*256 + 163], p7%2==0);
        compareSwap(arr[p7*256 + 164], arr[p7*256 + 165], arr_index[p7*256 + 164], arr_index[p7*256 + 165], p7%2==0);
        compareSwap(arr[p7*256 + 166], arr[p7*256 + 167], arr_index[p7*256 + 166], arr_index[p7*256 + 167], p7%2==0);
        compareSwap(arr[p7*256 + 168], arr[p7*256 + 169], arr_index[p7*256 + 168], arr_index[p7*256 + 169], p7%2==0);
        compareSwap(arr[p7*256 + 170], arr[p7*256 + 171], arr_index[p7*256 + 170], arr_index[p7*256 + 171], p7%2==0);
        compareSwap(arr[p7*256 + 172], arr[p7*256 + 173], arr_index[p7*256 + 172], arr_index[p7*256 + 173], p7%2==0);
        compareSwap(arr[p7*256 + 174], arr[p7*256 + 175], arr_index[p7*256 + 174], arr_index[p7*256 + 175], p7%2==0);
        compareSwap(arr[p7*256 + 176], arr[p7*256 + 177], arr_index[p7*256 + 176], arr_index[p7*256 + 177], p7%2==0);
        compareSwap(arr[p7*256 + 178], arr[p7*256 + 179], arr_index[p7*256 + 178], arr_index[p7*256 + 179], p7%2==0);
        compareSwap(arr[p7*256 + 180], arr[p7*256 + 181], arr_index[p7*256 + 180], arr_index[p7*256 + 181], p7%2==0);
        compareSwap(arr[p7*256 + 182], arr[p7*256 + 183], arr_index[p7*256 + 182], arr_index[p7*256 + 183], p7%2==0);
        compareSwap(arr[p7*256 + 184], arr[p7*256 + 185], arr_index[p7*256 + 184], arr_index[p7*256 + 185], p7%2==0);
        compareSwap(arr[p7*256 + 186], arr[p7*256 + 187], arr_index[p7*256 + 186], arr_index[p7*256 + 187], p7%2==0);
        compareSwap(arr[p7*256 + 188], arr[p7*256 + 189], arr_index[p7*256 + 188], arr_index[p7*256 + 189], p7%2==0);
        compareSwap(arr[p7*256 + 190], arr[p7*256 + 191], arr_index[p7*256 + 190], arr_index[p7*256 + 191], p7%2==0);

        compareSwap(arr[p7*256 + 192], arr[p7*256 + 193], arr_index[p7*256 + 192], arr_index[p7*256 + 193], p7%2==0);
        compareSwap(arr[p7*256 + 194], arr[p7*256 + 195], arr_index[p7*256 + 194], arr_index[p7*256 + 195], p7%2==0);
        compareSwap(arr[p7*256 + 196], arr[p7*256 + 197], arr_index[p7*256 + 196], arr_index[p7*256 + 197], p7%2==0);
        compareSwap(arr[p7*256 + 198], arr[p7*256 + 199], arr_index[p7*256 + 198], arr_index[p7*256 + 199], p7%2==0);
        compareSwap(arr[p7*256 + 200], arr[p7*256 + 201], arr_index[p7*256 + 200], arr_index[p7*256 + 201], p7%2==0);
        compareSwap(arr[p7*256 + 202], arr[p7*256 + 203], arr_index[p7*256 + 202], arr_index[p7*256 + 203], p7%2==0);
        compareSwap(arr[p7*256 + 204], arr[p7*256 + 205], arr_index[p7*256 + 204], arr_index[p7*256 + 205], p7%2==0);
        compareSwap(arr[p7*256 + 206], arr[p7*256 + 207], arr_index[p7*256 + 206], arr_index[p7*256 + 207], p7%2==0);
        compareSwap(arr[p7*256 + 208], arr[p7*256 + 209], arr_index[p7*256 + 208], arr_index[p7*256 + 209], p7%2==0);
        compareSwap(arr[p7*256 + 210], arr[p7*256 + 211], arr_index[p7*256 + 210], arr_index[p7*256 + 211], p7%2==0);
        compareSwap(arr[p7*256 + 212], arr[p7*256 + 213], arr_index[p7*256 + 212], arr_index[p7*256 + 213], p7%2==0);
        compareSwap(arr[p7*256 + 214], arr[p7*256 + 215], arr_index[p7*256 + 214], arr_index[p7*256 + 215], p7%2==0);
        compareSwap(arr[p7*256 + 216], arr[p7*256 + 217], arr_index[p7*256 + 216], arr_index[p7*256 + 217], p7%2==0);
        compareSwap(arr[p7*256 + 218], arr[p7*256 + 219], arr_index[p7*256 + 218], arr_index[p7*256 + 219], p7%2==0);
        compareSwap(arr[p7*256 + 220], arr[p7*256 + 221], arr_index[p7*256 + 220], arr_index[p7*256 + 221], p7%2==0);
        compareSwap(arr[p7*256 + 222], arr[p7*256 + 223], arr_index[p7*256 + 222], arr_index[p7*256 + 223], p7%2==0);
        compareSwap(arr[p7*256 + 224], arr[p7*256 + 225], arr_index[p7*256 + 224], arr_index[p7*256 + 225], p7%2==0);
        compareSwap(arr[p7*256 + 226], arr[p7*256 + 227], arr_index[p7*256 + 226], arr_index[p7*256 + 227], p7%2==0);
        compareSwap(arr[p7*256 + 228], arr[p7*256 + 229], arr_index[p7*256 + 228], arr_index[p7*256 + 229], p7%2==0);
        compareSwap(arr[p7*256 + 230], arr[p7*256 + 231], arr_index[p7*256 + 230], arr_index[p7*256 + 231], p7%2==0);
        compareSwap(arr[p7*256 + 232], arr[p7*256 + 233], arr_index[p7*256 + 232], arr_index[p7*256 + 233], p7%2==0);
        compareSwap(arr[p7*256 + 234], arr[p7*256 + 235], arr_index[p7*256 + 234], arr_index[p7*256 + 235], p7%2==0);
        compareSwap(arr[p7*256 + 236], arr[p7*256 + 237], arr_index[p7*256 + 236], arr_index[p7*256 + 237], p7%2==0);
        compareSwap(arr[p7*256 + 238], arr[p7*256 + 239], arr_index[p7*256 + 238], arr_index[p7*256 + 239], p7%2==0);
        compareSwap(arr[p7*256 + 240], arr[p7*256 + 241], arr_index[p7*256 + 240], arr_index[p7*256 + 241], p7%2==0);
        compareSwap(arr[p7*256 + 242], arr[p7*256 + 243], arr_index[p7*256 + 242], arr_index[p7*256 + 243], p7%2==0);
        compareSwap(arr[p7*256 + 244], arr[p7*256 + 245], arr_index[p7*256 + 244], arr_index[p7*256 + 245], p7%2==0);
        compareSwap(arr[p7*256 + 246], arr[p7*256 + 247], arr_index[p7*256 + 246], arr_index[p7*256 + 247], p7%2==0);
        compareSwap(arr[p7*256 + 248], arr[p7*256 + 249], arr_index[p7*256 + 248], arr_index[p7*256 + 249], p7%2==0);
        compareSwap(arr[p7*256 + 250], arr[p7*256 + 251], arr_index[p7*256 + 250], arr_index[p7*256 + 251], p7%2==0);
        compareSwap(arr[p7*256 + 252], arr[p7*256 + 253], arr_index[p7*256 + 252], arr_index[p7*256 + 253], p7%2==0);
        compareSwap(arr[p7*256 + 254], arr[p7*256 + 255], arr_index[p7*256 + 254], arr_index[p7*256 + 255], p7%2==0);
    }
}

