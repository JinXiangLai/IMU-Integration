//
// Created by hyj on 18-1-19.
//
#include <random>

#include "imu.h"
#include "utilities.h"

// euler2Rotation:  body frame to interitail frame
// IMU body系到惯性系下的旋转矩阵
Eigen::Matrix3d euler2Rotation( Eigen::Vector3d  eulerAngles)
{
    double roll = eulerAngles(0);
    double pitch = eulerAngles(1);
    double yaw = eulerAngles(2);

    double cr = cos(roll); double sr = sin(roll);
    double cp = cos(pitch); double sp = sin(pitch);
    double cy = cos(yaw); double sy = sin(yaw);

    Eigen::Matrix3d RIb;
    RIb << cy * cp, cy * sp * sr - sy * cr, sy * sr + cy * cr * sp, 
           sy * cp, cy * cr + sy * sr * sp, sp * sy * cr - cy * sr,
           -sp, cp * sr, cp * cr;
    return RIb;
}

// 欧拉角只涉及到旋转，这里是求
// 惯性系下的角速度到body系下的角速度，因为IMU测量值都是相对于它自身的坐标系测量的
Eigen::Matrix3d eulerRates2bodyRates(Eigen::Vector3d eulerAngles)
{
    double roll = eulerAngles(0);
    double pitch = eulerAngles(1);

    double cr = cos(roll); double sr = sin(roll);
    double cp = cos(pitch); double sp = sin(pitch);

    Eigen::Matrix3d R;
    R << 1, 0, -sp, 
         0, cr, sr * cp, 
         0, -sr, cr * cp;

    return R;
}

// IMU传感器的构造函数
IMU::IMU(Param p): param_(p)
{
    gyro_bias_ = Eigen::Vector3d::Zero();
    acc_bias_ = Eigen::Vector3d::Zero();
}

void IMU::addIMUnoise(MotionData& data)
{
    std::random_device rd;
    std::default_random_engine generator_(rd());
    // 0 均值的正态分布
    std::normal_distribution<double> noise(0.0, 1.0);

    Eigen::Vector3d noise_gyro(noise(generator_),noise(generator_),noise(generator_));
    Eigen::Matrix3d gyro_sqrt_cov = param_.gyro_noise_sigma * Eigen::Matrix3d::Identity();
    // 给IMU的角速度加上 高斯白噪声(和采样时间Δt成反比)，同时还要加上bias
    data.imu_gyro = data.imu_gyro + gyro_sqrt_cov * noise_gyro / sqrt( param_.imu_timestep ) + gyro_bias_;

    Eigen::Vector3d noise_acc(noise(generator_),noise(generator_),noise(generator_));
    Eigen::Matrix3d acc_sqrt_cov = param_.acc_noise_sigma * Eigen::Matrix3d::Identity();
    // 同样给IMU的加速度加上高斯白噪声和bias
    data.imu_acc = data.imu_acc + acc_sqrt_cov * noise_acc / sqrt( param_.imu_timestep ) + acc_bias_;

    // IMU的角速度和加速度对应的bias是随着时间进行随机游走的，
    // 这里利用噪声对其进行更新
    // gyro_bias update
    Eigen::Vector3d noise_gyro_bias(noise(generator_),noise(generator_),noise(generator_));
    gyro_bias_ += param_.gyro_bias_sigma * sqrt(param_.imu_timestep ) * noise_gyro_bias;
    data.imu_gyro_bias = gyro_bias_;

    // acc_bias update
    Eigen::Vector3d noise_acc_bias(noise(generator_),noise(generator_),noise(generator_));
    acc_bias_ += param_.acc_bias_sigma * sqrt(param_.imu_timestep ) * noise_acc_bias;
    data.imu_acc_bias = acc_bias_;

}

MotionData IMU::MotionModel(double t)
{

    MotionData data;
    // param
    float ellipse_x = 15;  // 椭圆短轴
    float ellipse_y = 20;  // 椭圆长轴
    float z = 1;           // z轴做sin运动
    float K1 = 10;          // z轴的正弦频率是x，y的k1倍
    float K = M_PI/ 10;    // 20 * K = 2pi 　　由于我们采取的是时间是20s, 系数K控制yaw正好旋转一圈，运动一周

    // translation 或者称为 Position
    // twb:  body frame in world frame
    // body系在世界系下 随时间 t 变化的位置
    Eigen::Vector3d position( ellipse_x * cos( K * t) + 5, ellipse_y * sin( K * t) + 5,  z * sin( K1 * K * t ) + 5);
    // body系在世界系下， 随时间t 变化的速度，直接用位置对时间的导数即可求到表达式
    Eigen::Vector3d dp(- K * ellipse_x * sin(K*t),  K * ellipse_y * cos(K*t), z*K1*K * cos(K1 * K * t));              // position导数　in world frame
    // body系在世界系下的加速度， 随时间t 变化， 使用位置求取二阶导数得到
    double K2 = K*K;
    Eigen::Vector3d ddp( -K2 * ellipse_x * cos(K*t),  -K2 * ellipse_y * sin(K*t), -z*K1*K1*K2 * sin(K1 * K * t));     // position二阶导数

    // Rotation
    double k_roll = 0.1;
    double k_pitch = 0.2;
    // 定义随时间变化的姿态， 这里使用欧拉角进行表示
    Eigen::Vector3d eulerAngles(k_roll * cos(t) , k_pitch * sin(t) , K*t );   // roll ~ [-0.2, 0.2], pitch ~ [-0.3, 0.3], yaw ~ [0,2pi]
    // 对姿态角求对时间的导数，得到其相对应的角速度
    Eigen::Vector3d eulerAnglesRates(-k_roll * sin(t) , k_pitch * cos(t) , K);      // euler angles 的导数

//    Eigen::Vector3d eulerAngles(0.0,0.0, K*t );   // roll ~ 0, pitch ~ 0, yaw ~ [0,2pi]
//    Eigen::Vector3d eulerAnglesRates(0.,0. , K);      // euler angles 的导数
    // eulerAngles是描述IMU body系自身的欧拉角，这里将其变换到世界系下的旋转矩阵
    Eigen::Matrix3d Rwb = euler2Rotation(eulerAngles);         // body frame to world frame
    // 三个欧拉角 r, p, y由于旋转有顺序，因此，其各轴下的角速度想要变换到Body系下，还要乘以矩阵 E_bw
    Eigen::Vector3d imu_gyro = eulerRates2bodyRates(eulerAngles) * eulerAnglesRates;   //  euler rates trans to body gyro
    
    // 东北天坐标系下的重力加速度
    Eigen::Vector3d gn (0,0,-9.81);                                   //  gravity in navigation frame(ENU)   ENU (0,0,-9.81)  NED(0,0,9,81)
    // IMU 坐标系下的加速度
    Eigen::Vector3d imu_acc = Rwb.transpose() * ( ddp -  gn );  //  Rbw * Rwn * gn = gs

    // IMU 的测量值只有这两个，
    // 陀螺仪测量的角速度
    data.imu_gyro = imu_gyro;
    // 加速度计测得的加速度
    data.imu_acc = imu_acc;

    // 这些量是通过对IMU 的gyro 和 acc进行离散积分得到的
    data.Rwb = Rwb;
    data.twb = position;
    data.imu_velocity = dp;
    data.timestamp = t;
    return data;

}

//读取生成的imu数据并用imu动力学模型对数据进行计算，最后保存imu积分以后的轨迹，
//用来验证数据以及模型的有效性。
// src表示生成的IMU data 即： imu_pose.txt or imu_pose_noise.txt
// dist表示根据IMU data生成的pose 即： imu_int_pose.txt or imu_int_pose_noise.txt (int 表示integration 积分)
void IMU::testImu(std::string src, std::string dist)
{
    // MotionData存储了IMU各个采样时刻的数据
    std::vector<MotionData>imudata;
    // 从外部文件中载入IMU数据
    LoadPose(src, imudata);

    std::ofstream save_points;
    save_points.open(dist);

    double dt = param_.imu_timestep;              // 定义采样步长Δt
    const double dt2 = dt * dt;
    Eigen::Vector3d Pwb = init_twb_;              // position :    from  imu measurements
    Eigen::Quaterniond Qwb(init_Rwb_);            // quaterniond:  from imu measurements
    Eigen::Vector3d Vw = init_velocity_;          // velocity  :   from imu measurements
    Eigen::Vector3d gw(0, 0, -9.81);                // ENU frame
    Eigen::Vector3d temp_a;
    Eigen::Vector3d theta;

    // 用于中值积分
    Eigen::Quaterniond Qwb0(init_Rwb_);

    // 用于预积分
    std::vector<Eigen::Vector3d> alpha_vector;
    std::vector<Eigen::Vector3d> belta_vector;
    std::vector<Eigen::Quaterniond> pre_q_vector;
    // 开始时刻，IMU没有进行运动，所以预积分应该都是0,因为Δt为0
    alpha_vector.push_back({ 0, 0, 0 });
    belta_vector.push_back({ 0, 0, 0 });
    pre_q_vector.push_back(Eigen::Quaterniond(1, 0, 0, 0));

    for (int i = 1; i < imudata.size(); ++i) {

        MotionData imupose0 = imudata[i - 1];
        MotionData imupose = imudata[i];

#define EULER  0
#define MID_POINT  0
#define PREINTEGRATION  1
#if EULER
        //delta_q = [1 , 1/2 * thetax , 1/2 * theta_y, 1/2 * theta_z]
        Eigen::Quaterniond dq;
        Eigen::Vector3d dtheta_half = imupose.imu_gyro * dt / 2.0;
        dq.w() = 1;
        dq.x() = dtheta_half.x();
        dq.y() = dtheta_half.y();
        dq.z() = dtheta_half.z();

        /// imu 动力学模型 欧拉积分
        // 世界系下的加速度
        Eigen::Vector3d acc_w = Qwb * (imupose.imu_acc) + gw;  // aw = Rwb * ( acc_body - acc_bias ) + gw
        // 更新IMU在世界系下的姿态
        Qwb = Qwb * dq;
        // IMU在世界系下的速度
        Vw = Vw + acc_w * dt;
        // IMU在世界系下的位置
        Pwb = Pwb + Vw * dt + 0.5 * dt2 * acc_w;

#elif MID_POINT
        /// 中值积分 ---> 需要先计算 Qk+1_w, 
        // 计算 Qk+1_w
        // 注意，imu_gyro 和 imu_acc数据
        // 这里的bias已经被融合到测量值里面，所以不需要减bias
        Eigen::Vector3d w_kplus1 = 0.5 * (imupose0.imu_gyro + imupose.imu_gyro);
        Eigen::Vector3d d_theta_half = w_kplus1 * dt * 0.5;
        Eigen::Quaterniond dq(1., d_theta_half.x(), d_theta_half.y(), d_theta_half.z());
        // 更新当前的姿态
        Qwb = Qwb0 * dq;

        Eigen::Vector3d acc_w = 0.5 * (Qwb0 * imupose0.imu_acc + Qwb * imupose.imu_acc) + gw;
        Vw = Vw + acc_w * dt;
        Pwb = Pwb + Vw * dt + 0.5 * acc_w * dt2;
        Qwb0 = Qwb;

#elif PREINTEGRATION
        // 预积分需要先把所有预积分算好后再一个个拿出来变到世界系下存储

        // IMU的角速度是绕其各轴的，所以参考系都是body系
        // 注意，要用临近的IMU测量进行中值计算，而不是和body0的测量数据
        Eigen::Vector3d w_i = 0.5 * (imupose0.imu_gyro + imupose.imu_gyro);
        // Eigen::Vector3d w_i = 0.5 * (imudata[0].imu_gyro + imupose.imu_gyro);

        // 算 q_b0_bi
        Eigen::Vector3d half_theta_i = 0.5 * w_i * dt;
        Eigen::Quaterniond delta_q(1., half_theta_i(0), half_theta_i(1), half_theta_i(2));
        Eigen::Quaterniond q_b0_bi = pre_q_vector[i - 1] * delta_q;
        pre_q_vector.push_back(q_b0_bi);

        // 算 a_b0_bi，注意，要用临近的IMU测量进行中值计算，而不是和body0的测量数据
        Eigen::Vector3d a_b0_bi = 0.5 * (pre_q_vector[i - 1] * imupose0.imu_acc + q_b0_bi * imupose.imu_acc);

        // 算预积分位移
        Eigen::Vector3d p_b0_bi = alpha_vector[i - 1] + belta_vector[i - 1] * dt + 0.5 * a_b0_bi * dt2;
        alpha_vector.push_back(p_b0_bi);

        // 计算速度
        Eigen::Vector3d v_b0_bi = belta_vector[i - 1] + a_b0_bi * dt;
        belta_vector.push_back(v_b0_bi);
    }
    // 把预积分数据转换到世界坐标系下
    const Eigen::Quaterniond q_wb0 = Eigen::Quaterniond(imudata[0].Rwb).normalized();
    const Eigen::Vector3d    t_wb0 = imudata[0].twb;
    const Eigen::Vector3d    v_wb0 = init_velocity_;
    for (int i = 1; i < imudata.size(); ++i) {
        MotionData imupose = imudata[i];
        const Eigen::Quaterniond Qwb = q_wb0 * pre_q_vector[i];
        // const Eigen::Vector3d    Pwb = t_wb0 + v_wb0 * (i*dt) + 0.5 * gw * dt2 + q_wb0 * alpha_vector[i];
        const Eigen::Vector3d    Pwb = t_wb0 + v_wb0 * (i*dt) + 0.5 * gw * (i * i * dt2) + q_wb0 * alpha_vector[i];
        save_points << imupose.timestamp << " "
            << Qwb.w() << " "
            << Qwb.x() << " "
            << Qwb.y() << " "
            << Qwb.z() << " "
            << Pwb(0) << " "
            << Pwb(1) << " "
            << Pwb(2) << " "
            << Qwb.w() << " "
            << Qwb.x() << " "
            << Qwb.y() << " "
            << Qwb.z() << " "
            << Pwb(0) << " "
            << Pwb(1) << " "
            << Pwb(2) << " "
            << std::endl;
    }
#endif

#if !PREINTEGRATION
        //　按着imu postion, imu quaternion , cam postion, cam quaternion 的格式存储，由于没有cam，所以imu存了两次
        // 把 P， V， Q存起来就好了
        save_points<<imupose.timestamp<<" "
                   <<Qwb.w()<<" "
                   <<Qwb.x()<<" "
                   <<Qwb.y()<<" "
                   <<Qwb.z()<<" "
                   <<Pwb(0)<<" "
                   <<Pwb(1)<<" "
                   <<Pwb(2)<<" "
                   <<Qwb.w()<<" "
                   <<Qwb.x()<<" "
                   <<Qwb.y()<<" "
                   <<Qwb.z()<<" "
                   <<Pwb(0)<<" "
                   <<Pwb(1)<<" "
                   <<Pwb(2)<<" "
                   <<std::endl;
    }
#endif

    std::cout<<"test　end"<<std::endl;

}
