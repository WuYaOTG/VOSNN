
#include "DSVS_Octree.h"
//_PL全局,k_为.h参数, hw是硬件要用的参数
static type_point data_set_xmin_PL;
static type_point data_set_xmax_PL;
static type_point data_set_ymin_PL;
static type_point data_set_ymax_PL;
static type_point data_set_zmin_PL;
static type_point data_set_zmax_PL;

static int split_x_array_size_PL;
static int split_y_array_size_PL;
static int split_z_array_size_PL;
static int total_voxel_size_PL;
static type_point voxel_split_unit_PL;


// 从 txt 读取点云到数组中
void read_points_from_txt(std::string file_name, My_PointXYZI* laserCloudInArray, int & point_size)
{
    std::ifstream file(file_name);
    if(!file)
    {
        std::cout << "!!! error! fails read pointcloud txt file from" << file_name << std::endl;
        point_size = 0;
    }
    else{
        std::string line;   // 读到的每一行数据
        std::string read_element;   // 每一行按空格拆分后的单个数据
        float line_data[6]; // 一行最多6个数据
        int line_i=0;   // 行数
        
        while(getline(file,line))   // 读取一行数据，读入line中
        {
            std::istringstream line_stream(line);   // 分割line中数据，转化成数据流
            int line_element_count = 0;
            while(line_stream >> read_element)  // 对于每一个数据，写入 line_data 中
            {
                line_data[line_element_count] = atof(read_element.data());
                line_element_count++;
            }

            {
                laserCloudInArray[line_i].x = line_data[0];
                laserCloudInArray[line_i].y = line_data[1];
                laserCloudInArray[line_i].z = line_data[2];
                // laserCloudInArray[line_i].intensity = line_data[3];s
            }

            line_i = line_i + 1;

        }

        point_size = line_i;

    }
}


void get_min_max(My_PointXYZI* data_set, int data_set_size, type_point_hw& data_set_xmin_hw, type_point_hw& data_set_ymin_hw, type_point_hw& data_set_zmin_hw)
{
	type_point x_max, x_min, y_max, y_min, z_max, z_min;

	x_max = y_max = z_max = -k_data_max_value_abs;
	x_min = y_min = z_min = k_data_max_value_abs;

    loop_cache_and_min_max:
	for (int i = 0; i < data_set_size; i++)
	{
		type_point data_x = data_set[i].x;
		type_point data_y = data_set[i].y;
		type_point data_z = data_set[i].z;

		//compare x
		if (x_max < data_x)
			x_max = data_x;
		if (x_min > data_x)
			x_min = data_x;

		//compare y
		if (y_max < data_y)
			y_max = data_y;
		if (y_min > data_y)
			y_min = data_y;

		//compare z
		if (z_max < data_z)
			z_max = data_z;
		if (z_min > data_z)
			z_min = data_z;
	}
	//store these min max in a structure for convenient usage
	//My_MaxMin data_set_max_min;
	data_set_xmin_PL = x_min;
	data_set_xmax_PL = x_max;
	data_set_ymin_PL = y_min;
	data_set_ymax_PL = y_max;
	data_set_zmin_PL = z_min;
	data_set_zmax_PL = z_max;

    data_set_xmin_hw = x_min;
    data_set_ymin_hw = y_min;
    data_set_zmin_hw = z_min;
}

void split_voxel_AABB(type_point split_unit, int& split_x_size_PL_hw, int& split_y_size_PL_hw, int& split_z_size_PL_hw, int& total_voxel_size_PL_hw)
{
    type_point split_unit_local = split_unit;
    voxel_split_unit_PL = split_unit;

    int x_split_size = ((data_set_xmax_PL - data_set_xmin_PL) / split_unit_local) + 1;
	int y_split_size = ((data_set_ymax_PL - data_set_ymin_PL) / split_unit_local) + 1;
	int z_split_size = ((data_set_zmax_PL - data_set_zmin_PL) / split_unit_local) + 1;

    if (x_split_size > k_axis_voxel_max)
        x_split_size = k_axis_voxel_max;
    if (y_split_size > k_axis_voxel_max)
        y_split_size = k_axis_voxel_max;
    if (z_split_size > k_axis_voxel_max)
        z_split_size = k_axis_voxel_max;

    split_x_array_size_PL = x_split_size;
	split_y_array_size_PL = y_split_size;
	split_z_array_size_PL = z_split_size;

    split_x_size_PL_hw = x_split_size;
	split_y_size_PL_hw = y_split_size;
	split_z_size_PL_hw = z_split_size;

    total_voxel_size_PL = x_split_size * y_split_size * z_split_size;
    total_voxel_size_PL_hw = total_voxel_size_PL;
}

void setup_hardware_PL(My_PointXYZI* KNN_reference_set, int reference_set_size,
    type_point_hw& data_set_xmin_hw, type_point_hw& data_set_ymin_hw, type_point_hw& data_set_zmin_hw, type_point_hw& voxel_split_unit_hw,
    int& split_x_size_PL_hw, int& split_y_size_PL_hw, int& split_z_size_PL_hw, int& total_voxel_size_PL_hw)
{
    get_min_max(KNN_reference_set, reference_set_size, data_set_xmin_hw, data_set_ymin_hw, data_set_zmin_hw);

    voxel_split_unit_hw = k_DSVS_split_unit;

    split_voxel_AABB(voxel_split_unit_hw, split_x_size_PL_hw, split_y_size_PL_hw, split_z_size_PL_hw, total_voxel_size_PL_hw);
}

//计算了每个点的hash值，由于是int，因此以1为单位的点会有相同hash值
void calculate_hash(My_PointXYZI* KNN_reference_set, int* data_set_hash, int KNN_reference_set_size)
{
    for (int i = 0; i < KNN_reference_set_size; i++)
	{
        //copy the point_data to data_x,y,z
        type_point data_x = KNN_reference_set[i].x;
        type_point data_y = KNN_reference_set[i].y;
        type_point data_z = KNN_reference_set[i].z;
        int x_split_array_size = split_x_array_size_PL;
        int y_split_array_size = split_y_array_size_PL;
        int z_split_array_size = split_z_array_size_PL;
        //default x,y,z index as the max index
        int x_index;
        int y_index;
        int z_index;

        if (data_x <= data_set_xmin_PL)
            x_index = 0;
        else
            x_index = (int)((data_x - data_set_xmin_PL) / voxel_split_unit_PL);
        if (x_index >= x_split_array_size)
            x_index = x_split_array_size - 1;

        if (data_y <= data_set_ymin_PL)
            y_index = 0;
        else
            y_index = (int)((data_y - data_set_ymin_PL) / voxel_split_unit_PL);
        if (y_index >= y_split_array_size)
            y_index = y_split_array_size - 1;

        if (data_z <= data_set_zmin_PL)
            z_index = 0;
        else
            z_index = (int)((data_z - data_set_zmin_PL) / voxel_split_unit_PL);
        if (z_index >= z_split_array_size)
            z_index = z_split_array_size - 1;

        //transform 3d index to a 1d index
        int data_hash = x_index * y_split_array_size * z_split_array_size + y_index * z_split_array_size + z_index;
        //因为数据集的hash值是从0开始的，因此需要减去1，因此也不会超过total_voxel_size_PL
        if (data_hash >= total_voxel_size_PL)
            data_hash = (total_voxel_size_PL - 1);
        if (data_hash < 0)
            data_hash = 0;

        data_set_hash[i] = data_hash;
    }
}

void DSVS_count_point(int* data_hash, int* DSVS_point_num, int KNN_reference_set_size)
{
    for (int i = 0; i < KNN_reference_set_size; i++)
	{
        int temp_hash = data_hash[i];
        DSVS_point_num[temp_hash] = DSVS_point_num[temp_hash] + 1;
    }
}

void DSVS_cal_first_index(int* DSVS_point_num, int* DSVS_first_index, int total_voxel_size)
{
    DSVS_first_index[0] = 0;
    for (int temp_hash = 1; temp_hash < total_voxel_size; temp_hash++)
	{
        DSVS_first_index[temp_hash] = DSVS_first_index[temp_hash - 1] + DSVS_point_num[temp_hash - 1];
    }
}

void DSVS_cal_flag(int* DSVS_first_index, int* DSVS_voxel_flag)
{
    for(int temp_hash = 0; temp_hash < total_voxel_size_PL; temp_hash++)
    {
        int voxel_point_size = DSVS_first_index[temp_hash + 1] - DSVS_first_index[temp_hash];
        if(voxel_point_size > k_Octree_threshold)
        {
            DSVS_voxel_flag[temp_hash] = k_div_flag[0];
        }
        else
        {
            DSVS_voxel_flag[temp_hash] = k_no_div_flag;
        }
    }
}

void DSVS_reorder(My_PointXYZI* set, int set_size, int* data_hash, int* DSVS_first_index, int* DSVS_point_num, int* count_point, 
    int* ordered_index, My_PointXYZI* ordered_point)
{
    for (int i = 0; i < set_size; i++)
    {
        int ordered_point_index;
        int temp_hash = data_hash[i];
        int start_point = DSVS_first_index[temp_hash];
        int end_point = start_point + DSVS_point_num[temp_hash];

        int count_point_temp = count_point[temp_hash];
        ordered_point_index = start_point + count_point_temp;
        count_point[temp_hash] = count_point_temp + 1;
        
        ordered_index[i] = ordered_point_index;
        ordered_point[ordered_point_index] = set[i];
    }
}

void Oc_depth0_init(int* Oc_voxel_depth, int* Oc_voxel_point_num, int* Oc_voxel_flag,
    int* Oc_voxel_first_index, int* DSVS_point_num, int* DSVS_first_index, int* DSVS_voxel_flag)
{
    for(int temp_hash = 0; temp_hash < total_voxel_size_PL; temp_hash++)
    {
        Oc_voxel_depth[temp_hash * k_max_encode] = 0;
        Oc_voxel_point_num[temp_hash * k_max_encode] = 0;
        Oc_voxel_flag[temp_hash * k_max_encode] = DSVS_voxel_flag[temp_hash];
        Oc_voxel_first_index[temp_hash * k_max_encode] = DSVS_first_index[temp_hash];
    }
}

void DSVS_build(My_PointXYZI* KNN_reference_set, int reference_set_size, int* DSVS_point_num, int* DSVS_first_index, 
    int* DSVS_voxel_flag, My_PointXYZI* DSVS_ordered_set, int* DSVS_ordered_index)
{
    int DSVS_hash[reference_set_size];
    int count_point_v[k_voxels_number_max] = {0};

    calculate_hash(KNN_reference_set, DSVS_hash, reference_set_size);
    DSVS_count_point(DSVS_hash, DSVS_point_num, reference_set_size);
    DSVS_cal_first_index(DSVS_point_num, DSVS_first_index, total_voxel_size_PL);
    DSVS_cal_flag(DSVS_first_index, DSVS_voxel_flag);
    DSVS_reorder(KNN_reference_set, reference_set_size, DSVS_hash, DSVS_first_index, DSVS_point_num, count_point_v, DSVS_ordered_index, DSVS_ordered_set);
}

void dep1_Oc_encode(My_PointXYZI* ordered_point, int* DSVS_voxel_flag, int* DSVS_first_index, int *DSVS_point_num, point_information* point_inform)
{
    for(int temp_hash = 0; temp_hash < total_voxel_size_PL; temp_hash++)
    {
        if(DSVS_voxel_flag[temp_hash] == k_div_flag[0])
        {
            int start_index = DSVS_first_index[temp_hash];
            // int end_index = start_index + Oc_point_num[temp_hash][0][0];
            
            int x_gain = split_y_array_size_PL * split_z_array_size_PL;
            int x_index = temp_hash / x_gain;
            int y_gain = split_z_array_size_PL;
            int y_index = (temp_hash - x_index * x_gain) / y_gain;
            int z_index = temp_hash - x_index * x_gain - y_index * y_gain;
            type_point xmin = data_set_xmin_PL + x_index * voxel_split_unit_PL;//type_point xmin = data_set_xmin_PL + x_index;
            type_point ymin = data_set_ymin_PL + y_index * voxel_split_unit_PL;//type_point ymin = data_set_ymin_PL + y_index;
            type_point zmin = data_set_zmin_PL + z_index * voxel_split_unit_PL;//type_point zmin = data_set_zmin_PL + z_index;

            for(int j = 0; j < DSVS_point_num[temp_hash]; j++)
            {
                int code = 0;
                int ordered_index_temp = start_index + j;
                My_PointXYZI point_temp = ordered_point[ordered_index_temp];

                code |= (point_temp.x - xmin >= k_mid_split_unit) << 2; // x 正为 1，负为 0
                code |= (point_temp.y - ymin >= k_mid_split_unit) << 1; // y 正为 1，负为 0
                code |= (point_temp.z - zmin >= k_mid_split_unit) << 0; // z 正为 1，负为 0

                // point_inform[ordered_index_temp].hash = temp_hash;
                point_inform[ordered_index_temp].encode = code;
                point_inform[ordered_index_temp].flag = k_div_flag[1];
                point_inform[ordered_index_temp].depth = 1;
            }
        }
        else
        {
            for(int j = 0; j < DSVS_point_num[temp_hash]; j++)
            {
                int ordered_index_temp = DSVS_first_index[temp_hash] + j;
                point_inform[ordered_index_temp].encode = 0;
                point_inform[ordered_index_temp].flag = k_no_div_flag;
                point_inform[ordered_index_temp].depth = 0;
            }
        }
    }
}

void dep1_Oc_point_num(int* Oc_voxel_point_num, int* Oc_voxel_depth, int* DSVS_voxel_flag, int* DSVS_first_index, 
    int *DSVS_point_num, point_information* point_inform)
{
    for(int temp_hash = 0; temp_hash < total_voxel_size_PL; temp_hash++)
    {
        if(DSVS_voxel_flag[temp_hash] == k_div_flag[0])
        {
            int start_index = DSVS_first_index[temp_hash];
            for(int j = 0; j < DSVS_point_num[temp_hash]; j++)
            {
                int temp_index = start_index + j;
                int temp_code = point_inform[temp_index].encode;
                Oc_voxel_point_num[temp_hash*k_max_encode + temp_code] = Oc_voxel_point_num[temp_hash*k_max_encode + temp_code] + 1;
                Oc_voxel_depth[temp_hash*k_max_encode + temp_code] = 1;
            }
        }
    }
}

//new version
// void dep1_Oc_point_num(point_information* point_inform, int* DSVS_hash, int set_size, int* dep1_point_num)
// {
//     for(int i = 0; i < set_size; i++)
//     {
//         point_information temp_point_info = point_inform[i];
//         int temp_hash = DSVS_hash[i];

//         int temp_code = temp_point_info.encode;
//         int addr = temp_hash * k_max_encode + temp_code;
//         dep1_point_num[addr]++;
//     }
// }

//new version
// void dep1_Oc_first_index(int* dep1_point_num, int* dep1_first_index)
// {
//     dep1_first_index[0] = 0;
//     for(int i = 1; i < total_voxel_size_hw * k_max_encode; i++)
//     {
//         int temp_hash = i / k_max_encode;
//         int temp_code = i % k_max_encode;
//         dep1_first_index[i] = dep1_first_index[i - 1] + dep1_point_num[i - 1];
//         // DEBUG_LOG("dep1_first_index[" << temp_hash << "][" << temp_code << "] = " << dep1_first_index[i]);
//     }
// }

void cal_Oc_first_index(int Depth, int* Oc_voxel_flag, int* Oc_voxel_first_index, int* Oc_voxel_point_num)
{// if depth = 1
    DEBUG_OPERATION(int div_num = 0;);
    for(int temp_hash = 0; temp_hash < total_voxel_size_PL; temp_hash++)
    {
        for(int last_code = 0; last_code < k_max_encode; last_code++)
        {
            if(Oc_voxel_flag[temp_hash*k_max_encode + last_code] == k_div_flag[Depth-1])
            {
                int start_index = Oc_voxel_first_index[temp_hash*k_max_encode + last_code];
                int start_code = (last_code<<3) | 0;
                Oc_voxel_first_index[temp_hash*k_max_encode + start_code] = start_index;
                for(int k = 1; k < 8; k++)
                {
                    int next_code = start_code + k;//(last_code<<3) | k
                    Oc_voxel_first_index[temp_hash*k_max_encode + next_code] = Oc_voxel_first_index[temp_hash*k_max_encode + next_code-1] + Oc_voxel_point_num[temp_hash*k_max_encode + next_code-1];
                    DEBUG_OPERATION(
                        if(Oc_voxel_point_num[temp_hash*k_max_encode + next_code-1] > k_Octree_threshold)
                        {
                            div_num++;
                        }
                    );
                }
            }
        }
    }
    DEBUG_OPERATION(
        if(div_num > 0)
        {
            DEBUG_LOG("Depth: " << Depth << "; div_num: " << div_num);
        }
    );
}

void cal_Oc_voxel_flag(int Depth, int* Oc_voxel_flag, int* Oc_voxel_depth, int* Oc_voxel_point_num)
{// if depth = 1
    for(int temp_hash = 0; temp_hash < total_voxel_size_PL; temp_hash++)
    {
        for(int temp_code = 0; temp_code < k_max_encode; temp_code++)//这里还有优化空间
        {
            if(Oc_voxel_depth[temp_hash*k_max_encode + temp_code] == Depth)
            {
                int temp_point_num = Oc_voxel_point_num[temp_hash*k_max_encode + temp_code];
                if(temp_point_num > k_Octree_threshold)
                {
                    Oc_voxel_flag[temp_hash*k_max_encode + temp_code] = k_div_flag[Depth];
                }
                else
                {
                    Oc_voxel_flag[temp_hash*k_max_encode + temp_code] = k_no_div_flag;//no div 是不是也要[]分层？
                }
            }
        }
    }
}

void cal_Oc_encode(int Depth, int* Oc_voxel_flag, int* Oc_voxel_first_index, int* Oc_voxel_point_num,
    My_PointXYZI* ordered_point, point_information* point_inform)
{
    type_point child_voxel_split_size  = voxel_split_unit_PL;
    for(int k = 1; k < Depth; k++)
    {
        child_voxel_split_size  = child_voxel_split_size / 2; // Attention!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(k == Depth-1)
        {
            DEBUG_LOG("Depth:" << Depth << "; " << "child_voxel_split_size: " << child_voxel_split_size);//DEBUG Attention!!!
        }
    }

    for(int temp_hash = 0; temp_hash < total_voxel_size_PL; temp_hash++)
    {
        for(int temp_code = 0; temp_code < k_max_encode; temp_code++)
        {
            if(Oc_voxel_flag[temp_hash*k_max_encode + temp_code] == k_div_flag[Depth-1])
            {
                int start_index = Oc_voxel_first_index[temp_hash*k_max_encode + temp_code];
                //int end_index = start_index + Oc_voxel[temp_hash][temp_code].point_num;

                int x_gain = split_y_array_size_PL * split_z_array_size_PL;
                int x_index = temp_hash / x_gain;
                int y_gain = split_z_array_size_PL;
                int y_index = (temp_hash - x_index * x_gain) / y_gain;
                int z_index = temp_hash - x_index * x_gain - y_index * y_gain;

                // 计算当前层级的累积偏移 Attention!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                type_point x_offset = 0.0f, y_offset = 0.0f, z_offset = 0.0f;
                for (int k = 1; k <= Depth; k++) {
                    type_point layer_size = voxel_split_unit_PL / (1 << (k - 1));
                    int layer_shift = 3 * (Depth - k);
                    int layer_code = (temp_code >> layer_shift) & 0x7;
                    x_offset += (layer_code & 0x4) ? layer_size : 0;
                    y_offset += (layer_code & 0x2) ? layer_size : 0;
                    z_offset += (layer_code & 0x1) ? layer_size : 0;
                }

                type_point xmin_parent = data_set_xmin_PL + x_index * voxel_split_unit_PL + x_offset;
                type_point ymin_parent = data_set_ymin_PL + y_index * voxel_split_unit_PL + y_offset;
                type_point zmin_parent = data_set_zmin_PL + z_index * voxel_split_unit_PL + z_offset;

                for(int k = 0; k < Oc_voxel_point_num[temp_hash*k_max_encode + temp_code]; k++)
                {
                    int temp_index = start_index + k;
                    My_PointXYZI temp_point = ordered_point[temp_index];

                    // 计算子体素编码
                    type_point x_rel = (temp_point.x - xmin_parent) / child_voxel_split_size;
                    type_point y_rel = (temp_point.y - ymin_parent) / child_voxel_split_size;
                    type_point z_rel = (temp_point.z - zmin_parent) / child_voxel_split_size;

                    int next_code = 0;
                    next_code |= (x_rel >= k_mid_split_unit) << 2;
                    next_code |= (y_rel >= k_mid_split_unit) << 1;
                    next_code |= (z_rel >= k_mid_split_unit) << 0;

                    int final_code = (temp_code << 3) | next_code;  // 合并编码
                    point_inform[temp_index].encode = final_code;
                    point_inform[temp_index].flag = k_div_flag[Depth];
                    point_inform[temp_index].depth = Depth;
                }
            }
        }
    }
}

void depN_cal_Oc_voxel_point(int Depth, int* Oc_voxel_point_num, int* next_dep_point_num, int* Oc_voxel_first_index, 
    int* Oc_voxel_depth, int* Oc_voxel_flag, point_information* point_inform)
{
    for(int temp_hash = 0; temp_hash < total_voxel_size_PL; temp_hash++)
    {
        for(int temp_code = 0; temp_code < k_max_encode; temp_code++)
        {
            if(Oc_voxel_flag[temp_hash*k_max_encode + temp_code] == k_div_flag[Depth-1])
            {
                int start_index = Oc_voxel_first_index[temp_hash*k_max_encode + temp_code];
                for(int k = 0; k < Oc_voxel_point_num[temp_hash*k_max_encode + temp_code]; k++)
                {
                    int temp_index = start_index + k;
                    int temp_final_code = point_inform[temp_index].encode;
                    next_dep_point_num[temp_hash*k_max_encode + temp_final_code] = next_dep_point_num[temp_hash*k_max_encode + temp_final_code] + 1;
                    Oc_voxel_depth[temp_hash*k_max_encode + temp_final_code] = Depth;
                }
            }
        }
    }
}

void depN_Octree_reorder(int Depth, My_PointXYZI* set, int* set_index, int set_size, int* data_hash, My_PointXYZI* depthN_set, int* depthN_set_index,
    int* Oc_voxel_first_index, point_information* point_inform, point_information* depN_point_inform, int* count_temp)
{//if depth = 2
    for(int i = 0; i < set_size; i++)
    {
        int temp_hash = data_hash[i];
        int temp_flag = point_inform[i].flag;

        if(temp_flag == k_div_flag[Depth])
        {
            int temp_code = point_inform[i].encode;//int temp_code = point_inform[temp_index].encode;
            int start_index = Oc_voxel_first_index[temp_hash*k_max_encode + temp_code];

            int temp_ordered_index = start_index + count_temp[temp_hash*k_max_encode + temp_code];
            count_temp[temp_hash*k_max_encode + temp_code] = count_temp[temp_hash*k_max_encode + temp_code] + 1;

            depthN_set_index[i] = temp_ordered_index;
            depthN_set[temp_ordered_index] = set[i];
            depN_point_inform[temp_ordered_index].encode = point_inform[i].encode;
        }
        else
        {
            depthN_set_index[i] = set_index[i];
            depthN_set[i] = set[i];
            depN_point_inform[i].encode = point_inform[i].encode;
        }
    }
}

void Octree_Init(My_PointXYZI* DSVS_ordered_set, int* DSVS_ordered_index, int set_size, int* DSVS_point_num, int* DSVS_first_index,
    int* DSVS_voxel_flag, int* Oc_voxel_depth, int* Oc_voxel_point_num, int* Oc_voxel_flag, int* Oc_voxel_first_index,
    point_information* point_inform, point_information* dep1_point_inform, My_PointXYZI* dep1_set, int* dep1_set_index)
{
    int dep0_hash[set_size];
    int depN1_count[k_voxels_number_max * k_max_encode] = {0};

    Oc_depth0_init(Oc_voxel_depth, Oc_voxel_point_num, Oc_voxel_flag, Oc_voxel_first_index, DSVS_point_num, DSVS_first_index, DSVS_voxel_flag);
    calculate_hash(DSVS_ordered_set, dep0_hash, set_size);
    dep1_Oc_encode(DSVS_ordered_set, DSVS_voxel_flag, DSVS_first_index, DSVS_point_num, point_inform);
    dep1_Oc_point_num(Oc_voxel_point_num, Oc_voxel_depth, DSVS_voxel_flag, DSVS_first_index, DSVS_point_num, point_inform);
    cal_Oc_first_index(1, Oc_voxel_flag, Oc_voxel_first_index, Oc_voxel_point_num);
    cal_Oc_voxel_flag(1, Oc_voxel_flag, Oc_voxel_depth, Oc_voxel_point_num);
    depN_Octree_reorder(1, DSVS_ordered_set, DSVS_ordered_index, set_size, dep0_hash, dep1_set, dep1_set_index, Oc_voxel_first_index, point_inform, dep1_point_inform, depN1_count);
}

void Octree_build(int Depth, My_PointXYZI* dep1_set, int* dep1_set_index, int set_size, int* Oc_voxel_depth, int* Oc_voxel_point_num,
    int* Oc_voxel_flag, int* Oc_voxel_first_index, point_information* dep1_point_inform, point_information* dep2_point_inform,
    int* dep2_Oc_point_num, My_PointXYZI* dep2_set, int* dep2_set_index)
{
    int dep1_hash[set_size];
    int depN2_count[k_voxels_number_max * k_max_encode] = {0};

    calculate_hash(dep1_set, dep1_hash, set_size);
    cal_Oc_encode(Depth, Oc_voxel_flag, Oc_voxel_first_index, Oc_voxel_point_num, dep1_set, dep1_point_inform);
    depN_cal_Oc_voxel_point(Depth, Oc_voxel_point_num, dep2_Oc_point_num, Oc_voxel_first_index, Oc_voxel_depth, Oc_voxel_flag, dep1_point_inform);
    cal_Oc_first_index(Depth, Oc_voxel_flag, Oc_voxel_first_index, dep2_Oc_point_num);
    cal_Oc_voxel_flag(Depth, Oc_voxel_flag, Oc_voxel_depth, dep2_Oc_point_num);
    depN_Octree_reorder(Depth, dep1_set, dep1_set_index, set_size, dep1_hash, dep2_set, dep2_set_index, Oc_voxel_first_index, dep1_point_inform, dep2_point_inform, depN2_count);
}

void Octree_test(My_PointXYZI* DSVS_ordered_set, int* DSVS_ordered_index, int set_size, int* DSVS_point_num, int* DSVS_first_index,
    int* DSVS_voxel_flag, int* Oc_voxel_depth, int* Oc_depN_point_num, int* Oc_voxel_flag, int* Oc_voxel_first_index, 
    point_information* point_inform, point_information* depN_point_inform, My_PointXYZI* depN_set, int* depN_set_index)
{
    int Depth;
    point_information dep1_point_inform[set_size];
    point_information dep2_point_inform[set_size];
    My_PointXYZI dep1_set[set_size];
    My_PointXYZI dep2_set[set_size];
    int dep1_set_index[set_size];
    int dep2_set_index[set_size];
    int Oc_dep1_point_num[k_voxels_number_max * k_max_encode] = {0};
    int Oc_dep2_point_num[k_voxels_number_max * k_max_encode] = {0};

    Octree_Init(DSVS_ordered_set, DSVS_ordered_index, set_size, DSVS_point_num, DSVS_first_index, DSVS_voxel_flag, Oc_voxel_depth,
        Oc_dep1_point_num, Oc_voxel_flag, Oc_voxel_first_index, point_inform, dep1_point_inform, dep1_set, dep1_set_index);

    Octree_build(2, dep1_set, dep1_set_index, set_size, Oc_voxel_depth, Oc_dep1_point_num, Oc_voxel_flag, Oc_voxel_first_index,
        dep1_point_inform, depN_point_inform, Oc_depN_point_num, depN_set, depN_set_index);

    // Octree_build(3, dep2_set, dep2_set_index, set_size, Oc_voxel_depth, Oc_dep2_point_num, Oc_voxel_flag, Oc_voxel_first_index,
    //     dep2_point_inform, depN_point_inform, Oc_depN_point_num, depN_set, depN_set_index);
}

void cal_query_point_inform_new(My_PointXYZI* query_set, int q_set_size, int* q_hash, int* DSVS_voxel_flag, int* Oc_voxel_flag, 
    int* DSVS_first_index, int* Oc_voxel_first_index, int* q_flag, voxel_int* q_set_start_index)
{
    type_point child_voxel_split_size = voxel_split_unit_PL / 2; // Attention!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for(int i = 0; i < q_set_size; i++)
    {
        int temp_hash = q_hash[i];
        q_flag[i] = k_no_div_flag;

        int x_gain = split_y_array_size_PL * split_z_array_size_PL;
        int x_index = temp_hash / x_gain;
        int y_gain = split_z_array_size_PL;
        int y_index = (temp_hash - x_index * x_gain) / y_gain;
        int z_index = temp_hash - x_index * x_gain - y_index * y_gain;

        type_point xmin = data_set_xmin_PL + x_index * voxel_split_unit_PL;
        type_point ymin = data_set_ymin_PL + y_index * voxel_split_unit_PL;
        type_point zmin = data_set_zmin_PL + z_index * voxel_split_unit_PL;

        int code = 0;
        code |= (query_set[i].x - xmin >= k_mid_split_unit) << 2;
        code |= (query_set[i].y - ymin >= k_mid_split_unit) << 1;
        code |= (query_set[i].z - zmin >= k_mid_split_unit) << 0;

        if(DSVS_voxel_flag[temp_hash] == k_div_flag[0])
        {
            q_flag[i] = k_div_flag[0];
            q_set_start_index[i] = Oc_voxel_first_index[temp_hash*k_max_encode + code];
            if(Oc_voxel_flag[temp_hash*k_max_encode + code] == k_div_flag[1])
            {
                q_flag[i] = k_div_flag[1];

                type_point xmin_child = xmin + (code & 0x4 ? child_voxel_split_size : 0);
                type_point ymin_child = ymin + (code & 0x2 ? child_voxel_split_size : 0);
                type_point zmin_child = zmin + (code & 0x1 ? child_voxel_split_size : 0);

                type_point x_rel = (query_set[i].x - xmin_child) / child_voxel_split_size;
                type_point y_rel = (query_set[i].y - ymin_child) / child_voxel_split_size;
                type_point z_rel = (query_set[i].z - zmin_child) / child_voxel_split_size;

                int sub_code = 0;
                sub_code |= (x_rel >= k_mid_split_unit) << 2;
                sub_code |= (y_rel >= k_mid_split_unit) << 1;
                sub_code |= (z_rel >= k_mid_split_unit) << 0;

                int final_code = (code << 3) | sub_code;
                q_set_start_index[i] = Oc_voxel_first_index[temp_hash*k_max_encode + final_code];
            }
        }
        else
        {
            q_set_start_index[i] = DSVS_first_index[temp_hash];
        }

    }
}

void cal_query_info_v2(int q_set_size, int* q_hash, int* DSVS_voxel_flag,
    int* DSVS_first_index, int* q_flag, voxel_int* q_set_start_index)
{
    for(int i = 0; i < q_set_size; i++)
    {
        int temp_hash = q_hash[i];
        q_flag[i] = k_no_div_flag;

        if(DSVS_voxel_flag[temp_hash] == k_div_flag[0])
        {
            q_flag[i] = k_div_flag[0];
            q_set_start_index[i] = DSVS_first_index[temp_hash];
        }
        else
        {
            q_flag[i] = k_no_div_flag;
            q_set_start_index[i] = DSVS_first_index[temp_hash];
        }
    }
}

void query_init_v2(My_PointXYZI* query_set, int q_set_size, int* DSVS_voxel_flag, int* DSVS_first_index, 
    int* q_flag, voxel_int* q_set_start_index)
{
    int q_hash[q_set_size];

    calculate_hash(query_set, q_hash, q_set_size);
    cal_query_info_v2(q_set_size, q_hash, DSVS_voxel_flag, DSVS_first_index, q_flag, q_set_start_index);
}

void query_init(My_PointXYZI* query_set, int q_set_size, int* DSVS_voxel_flag, int* DSVS_first_index, 
    int* Oc_voxel_first_index, int* Oc_voxel_flag, int* q_flag, voxel_int* q_set_start_index)
{
    int q_hash[q_set_size];

    calculate_hash(query_set, q_hash, q_set_size);
    cal_query_point_inform_new(query_set, q_set_size, q_hash, DSVS_voxel_flag, Oc_voxel_flag, DSVS_first_index, Oc_voxel_first_index, q_flag, q_set_start_index);
}

// 暴力搜索
void brute_force_search(My_PointXYZI* query_set, int query_set_size, My_PointXYZI* reference_set, int reference_set_size, My_PointXYZI* bf_KNN_query_result)
{
    TIMER_INIT(5);
    TIMER_START(0);

    for(int query_i = 0; query_i < query_set_size; query_i++)
    {
        My_PointXYZI query_point = query_set[query_i];
        My_PointXYZI this_query_KNN_result;
        this_query_KNN_result.x = 0; this_query_KNN_result.y = 0;  this_query_KNN_result.z = 0;  //this_query_KNN_result.intensity = 0; 

        float current_min_dis = 10000;

        for(int i = 0; i < reference_set_size; i++)
        {
            float this_dis = sqrt( 
                (query_point.x- reference_set[i].x) * (query_point.x- reference_set[i].x) + 
                (query_point.y- reference_set[i].y) * (query_point.y- reference_set[i].y) + 
                (query_point.z- reference_set[i].z) * (query_point.z- reference_set[i].z) 
            );

            if(this_dis < current_min_dis)
            {
                this_query_KNN_result.x = reference_set[i].x;
                this_query_KNN_result.y = reference_set[i].y;
                this_query_KNN_result.z = reference_set[i].z;
                // this_query_KNN_result.intensity = reference_set[i].intensity;
                current_min_dis = this_dis;
            }
        }

        bf_KNN_query_result[query_i] = this_query_KNN_result;
    }
    int count = 0;

    TIMER_STOP_ID(0);

    DEBUG_TIME("Finished brute force search with query points " << query_set_size << " with " << TIMER_REPORT_MS(0) << " ms !" );
}

// 比较结果
void compare_result(My_PointXYZI* query_set, int query_set_size, My_PointXYZI* new_query_result, My_PointXYZI* gt_query_result, int & error_count)
{
    int count = 0;
    for(int i = 0; i < query_set_size; i++)
    {
        if(fabs(gt_query_result[i].x - new_query_result[i].x) > 0.1 || 
            fabs(gt_query_result[i].y - new_query_result[i].y) > 0.1 || 
            fabs(gt_query_result[i].z - new_query_result[i].z) > 0.1 )
        {
            count = count + 1;

            float dis_gt = sqrt( 
                (query_set[i].x- gt_query_result[i].x) * (query_set[i].x- gt_query_result[i].x) + 
                (query_set[i].y- gt_query_result[i].y) * (query_set[i].y- gt_query_result[i].y) + 
                (query_set[i].z- gt_query_result[i].z) * (query_set[i].z- gt_query_result[i].z) 
            );
            float dis_new = sqrt( 
                (query_set[i].x- new_query_result[i].x) * (query_set[i].x- new_query_result[i].x) + 
                (query_set[i].y- new_query_result[i].y) * (query_set[i].y- new_query_result[i].y) + 
                (query_set[i].z- new_query_result[i].z) * (query_set[i].z- new_query_result[i].z) 
            );

            if(dis_gt < 1)//fabs(dis_gt - dis_new) > 0.1
            {
                error_count ++;
                // int temp_hash = q_point_inform[i].hash;
                // int temp_encode = q_point_inform[i].encode;
                // int temp_flag = q_point_inform[i].flag;
                // DEBUG_INFO("Error, query flag: " << temp_flag << ", hash: " << temp_hash << ", encode: " << temp_encode);
            }
        }
    }
    DEBUG_INFO("compare_result count: " << count);
}
