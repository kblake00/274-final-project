#include "mbed.h"
#include "rtos.h"
#include "EthernetInterface.h"
#include "ExperimentServer.h"
#include "QEI.h"
#include "BezierCurve.h"
#include "MotorShield.h" 
#include "HardwareSetup.h"
#include "Matrix.h"
#include "MatrixMath.h"

#define BEZIER_ORDER_FOOT 7
#define NUM_INPUTS 22 // changed from (14 + 2*(BEZIER_ORDER_FOOT+1)) b/c no longer using bezier
#define NUM_OUTPUTS 45

#define PULSE_TO_RAD (2.0f*3.14159f / 1200.0f)

// Initializations
Serial pc(USBTX, USBRX);    // USB Serial Terminal
ExperimentServer server;    // Object that lets us communicate with MATLAB
Timer t;                    // Timer to measure elapsed time of experiment

QEI encoderA(PE_9,PE_11, NC, 1200, QEI::X4_ENCODING);  // MOTOR A ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderB(PA_5, PB_3, NC, 1200, QEI::X4_ENCODING);  // MOTOR B ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderC(PC_6, PC_7, NC, 1200, QEI::X4_ENCODING);  // MOTOR C ENCODER (no index, 1200 counts/rev, Quadrature encoding)
QEI encoderD(PD_12, PD_13, NC, 1200, QEI::X4_ENCODING);// MOTOR D ENCODER (no index, 1200 counts/rev, Quadrature encoding)

MotorShield motorShield(24000); //initialize the motor shield with a period of 12000 ticks or ~20kHZ
Ticker currentLoop;

Matrix MassMatrix(2,2);
Matrix JacobianL(2,2);
Matrix JacobianLT(2,2);
Matrix JacobianR(2,2);
Matrix JacobianRT(2,2);
Matrix tauL(2,1);
Matrix tauR(2,1);

// q1 is hip joint on non-inverted (right) leg -> motor A
// q2 is knee joint on non-inverted (right) leg -> motor B
// q3 is hip joint on inverted (left) leg -> motor C
// q4 is knee joint on inverted (left) leg -> motor D
// non-inverted leg is the leg closest to you if viewing the robot walking forward corresponds to the robot moving to the right

// Variables for q1
float current1;
float current_des1 = 0;
float prev_current_des1 = 0;
float current_int1 = 0;
float angle1;
float velocity1;
float duty_cycle1;
float angle1_init;

// Variables for q2
float current2;
float current_des2 = 0;
float prev_current_des2 = 0;
float current_int2 = 0;
float angle2;
float velocity2;
float duty_cycle2;
float angle2_init;

// Variables for q3
float current3;
float current_des3 = 0;
float prev_current_des3 = 0;
float current_int3 = 0;
float angle3;
float velocity3;
float duty_cycle3;
float angle3_init;

// Variables for q4
float current4;
float current_des4 = 0;
float prev_current_des4 = 0;
float current_int4 = 0;
float angle4;
float velocity4;
float duty_cycle4;
float angle4_init;

// Fixed kinematic parameters
const float l_OA=.011; 
const float l_OB=.042; 
const float l_AC=.096; 
const float l_DE=.091;
const float m1 =.0393 + .2;
const float m2 =.0368; 
const float m3 = .00783;
const float m4 = .0155;
const float I1 = 0.0000251;  //25.1 * 10^-6;
const float I2 = 0.0000535;  //53.5 * 10^-6;
const float I3 = 0.00000925; //9.25 * 10^-6;
const float I4 = 0.0000222;  //22.176 * 10^-6;
const float l_O_m1=0.032;
const float l_B_m2=0.0344; 
const float l_A_m3=0.0622;
const float l_C_m4=0.0610;
const float N = 18.75;
const float Ir = 0.0035/pow(N,2);

// Timing parameters
float current_control_period_us = 200.0f;     // 5kHz current control loop
float impedance_control_period_us = 1000.0f;  // 1kHz impedance control loop
float start_period, traj_period, end_period;

// Control parameters
float current_Kp = 4.0f;         
float current_Ki = 0.4f;           
float current_int_max = 3.0f;       
float duty_max;      
float K_xx;
float K_yy;
float K_xy;
float D_xx;
float D_xy;
float D_yy;

// trajectory parameters for left leg
float omega; 
float x_0;
float y_0;
float r;
float offset;

//ellipsoid parameters
float a;
float b;
float theta;

// Model parameters
float supply_voltage = 12;     // motor supply voltage
float R = 2.0f;                // motor resistance
float k_t = 0.18f;             // motor torque constant
float nu = 0.0005;             // motor viscous friction

// Current control interrupt function
void CurrentLoop()
{
    // This loop sets the motor voltage commands using PI current controllers with feedforward terms.
    
    //use the motor shield as follows:
    //motorShield.motorAWrite(DUTY CYCLE, DIRECTION), DIRECTION = 0 is forward, DIRECTION =1 is backwards.
    
    // q1 motor current control
    current1 = -(((float(motorShield.readCurrentA())/65536.0f)*30.0f)-15.0f);           // measure current
    velocity1 = encoderA.getVelocity() * PULSE_TO_RAD;                                  // measure velocity        
    float err_c1 = current_des1 - current1;                                             // current errror
    current_int1 += err_c1;                                                             // integrate error
    current_int1 = fmaxf( fminf(current_int1, current_int_max), -current_int_max);      // anti-windup
    float ff1 = R*current_des1 + k_t*velocity1;                                         // feedforward terms
    duty_cycle1 = (ff1 + current_Kp*err_c1 + current_Ki*current_int1)/supply_voltage;   // PI current controller
    
    float absDuty1 = abs(duty_cycle1);
    if (absDuty1 > duty_max) {
        duty_cycle1 *= duty_max / absDuty1;
        absDuty1 = duty_max;
    }    
    if (duty_cycle1 < 0) { // backwards
        motorShield.motorAWrite(absDuty1, 1);
    } else { // forwards
        motorShield.motorAWrite(absDuty1, 0);
    }             
    prev_current_des1 = current_des1; 
    
    // q2 motor current control
    current2     = -(((float(motorShield.readCurrentB())/65536.0f)*30.0f)-15.0f);       // measure current
    velocity2 = encoderB.getVelocity() * PULSE_TO_RAD;                                  // measure velocity  
    float err_c2 = current_des2 - current2;                                             // current error
    current_int2 += err_c2;                                                             // integrate error
    current_int2 = fmaxf( fminf(current_int2, current_int_max), -current_int_max);      // anti-windup   
    float ff2 = R*current_des2 + k_t*velocity2;                                         // feedforward terms
    duty_cycle2 = (ff2 + current_Kp*err_c2 + current_Ki*current_int2)/supply_voltage;   // PI current controller
    
    float absDuty2 = abs(duty_cycle2);
    if (absDuty2 > duty_max) {
        duty_cycle2 *= duty_max / absDuty2;
        absDuty2 = duty_max;
    }    
    if (duty_cycle2 < 0) { // backwards
        motorShield.motorBWrite(absDuty2, 1);
    } else { // forwards
        motorShield.motorBWrite(absDuty2, 0);
    }             
    prev_current_des2 = current_des2; 
    
    // q3 motor current control
    current3 = -(((float(motorShield.readCurrentC())/65536.0f)*30.0f)-15.0f);           // measure current
    velocity3 = encoderC.getVelocity() * PULSE_TO_RAD;                                  // measure velocity        
    float err_c3 = current_des3 - current3;                                             // current errror
    current_int3 += err_c3;                                                             // integrate error
    current_int3 = fmaxf( fminf(current_int3, current_int_max), -current_int_max);      // anti-windup
    float ff3 = R*current_des3 + k_t*velocity3;                                         // feedforward terms
    duty_cycle3 = (ff3 + current_Kp*err_c3 + current_Ki*current_int3)/supply_voltage;   // PI current controller
    
    float absDuty3 = abs(duty_cycle3);
    if (absDuty3 > duty_max) {
        duty_cycle3 *= duty_max / absDuty3;
        absDuty3 = duty_max;
    }    
    if (duty_cycle3 < 0) { // backwards
        motorShield.motorCWrite(absDuty3, 1);
    } else { // forwards
        motorShield.motorCWrite(absDuty3, 0);
    }             
    prev_current_des3 = current_des3;

    // q4 motor current control
    current4 = -(((float(motorShield.readCurrentD())/65536.0f)*30.0f)-15.0f);           // measure current
    velocity4 = encoderD.getVelocity() * PULSE_TO_RAD;                                  // measure velocity        
    float err_c4 = current_des4 - current4;                                             // current errror
    current_int4 += err_c4;                                                             // integrate error
    current_int4 = fmaxf( fminf(current_int4, current_int_max), -current_int_max);      // anti-windup
    float ff4 = R*current_des4 + k_t*velocity4;                                         // feedforward terms
    duty_cycle4 = (ff4 + current_Kp*err_c4 + current_Ki*current_int4)/supply_voltage;   // PI current controller
    
    float absDuty4 = abs(duty_cycle4);
    if (absDuty4 > duty_max) {
        duty_cycle4 *= duty_max / absDuty4;
        absDuty4 = duty_max;
    }    
    if (duty_cycle4 < 0) { // backwards
        motorShield.motorDWrite(absDuty4, 1);
    } else { // forwards
        motorShield.motorDWrite(absDuty4, 0);
    }             
    prev_current_des4 = current_des4;
}

int main (void)
{
    
    // Link the terminal with our server and start it up
    server.attachTerminal(pc);
    server.init();
    
    // Continually get input from MATLAB and run experiments
    float input_params[NUM_INPUTS];
    pc.printf("%f",input_params[0]);
    
    while(1) {
        
        // If there are new inputs, this code will run
        if (server.getParams(input_params,NUM_INPUTS)) {
           
                        
            // Get inputs from MATLAB          
            start_period                = input_params[0];    // First buffer time, before trajectory
            traj_period                 = input_params[1];    // Trajectory time/length
            end_period                  = input_params[2];    // Second buffer time, after trajectory
    
            angle1_init                 = input_params[3];    // Initial angle for q1 (rad)
            angle2_init                 = input_params[4];    // Initial angle for q2 (rad)
            angle3_init                 = input_params[5];    // Initial angle for q3 (rad)
            angle4_init                 = input_params[6];    // Initial angle for q4 (rad)

            K_xx                        = input_params[7];    // Foot stiffness N/m
            K_yy                        = input_params[8];    // Foot stiffness N/m
            K_xy                        = input_params[9];    // Foot stiffness N/m
            D_xx                        = input_params[10];    // Foot damping N/(m/s)
            D_yy                        = input_params[11];    // Foot damping N/(m/s)
            D_xy                        = input_params[12];   // Foot damping N/(m/s)
            duty_max                    = input_params[13];   // Maximum duty factor

            omega                     = input_params[14]; // trajectory parameters
            x_0                       = input_params[15];
            y_0                       = input_params[16];
            r                         = input_params[17];
            offset                      = input_params[18];

            a                           = input_params[19]; // ellipsoid parameters
            b                           = input_params[20];
            theta                       = input_params[21];
            
            // Attach current loop interrupt
            currentLoop.attach_us(CurrentLoop,current_control_period_us);
                        
            // Setup experiment
            t.reset();
            t.start();
            encoderA.reset();
            encoderB.reset();
            encoderC.reset();
            encoderD.reset();

            motorShield.motorAWrite(0, 0); //turn motor A off
            motorShield.motorBWrite(0, 0); //turn motor B off
            motorShield.motorCWrite(0, 0); //turn motor C off
            motorShield.motorDWrite(0, 0); //turn motor D off
                         
            // Run experiment
            while( t.read() < start_period + traj_period) {
                 
                // Read encoders to get motor states
                angle1 = encoderA.getPulses() *PULSE_TO_RAD + angle1_init;       
                velocity1 = encoderA.getVelocity() * PULSE_TO_RAD;
                 
                angle2 = encoderB.getPulses() * PULSE_TO_RAD + angle2_init;       
                velocity2 = encoderB.getVelocity() * PULSE_TO_RAD;   

                angle3 = encoderC.getPulses() *PULSE_TO_RAD + angle3_init;
                velocity3 = encoderC.getPulses() *PULSE_TO_RAD;

                angle4 = encoderD.getPulses() *PULSE_TO_RAD + angle4_init;
                velocity4 = encoderD.getPulses() *PULSE_TO_RAD;        

                const float th1 = angle1;
                const float th2 = angle2;
                const float th3 = angle3;
                const float th4 = angle4;
                const float dth1 = velocity1;
                const float dth2 = velocity2;
                const float dth3 = velocity3;
                const float dth4 = velocity4;
 
                // Calculate the Jacobian
                float Jx_th1 = l_AC*cos(th1 + th2) + l_DE*cos(th1) + l_OB*cos(th1);
                float Jx_th2 = l_AC*cos(th1 + th2);
                float Jy_th1 = l_AC*sin(th1 + th2) + l_DE*sin(th1) + l_OB*sin(th1);
                float Jy_th2 = l_AC*sin(th1 + th2);

                float Jx_th3 = l_AC*cos(th3 + th4) + l_DE*cos(th3) + l_OB*cos(th3);
                float Jx_th4 = l_AC*cos(th3 + th4);
                float Jy_th3 = l_AC*sin(th3 + th4) + l_DE*sin(th3) + l_OB*sin(th3);
                float Jy_th4 = l_AC*sin(th3 + th4);

                // Calculate the forward kinematics (position and velocity)
                float xFootR = l_AC*sin(th1 + th2) + l_DE*sin(th1) + l_OB*sin(th1);
                float yFootR = - l_AC*cos(th1 + th2) - l_DE*cos(th1) - l_OB*cos(th1);
                float dxFootR = Jx_th1*dth1 + Jx_th2*dth2;
                float dyFootR = Jy_th1*dth1 + Jy_th2*dth2;

                float xFootL = l_AC*sin(th3 + th4) + l_DE*sin(th3) + l_OB*sin(th3);
                float yFootL = - l_AC*cos(th3 + th4) - l_DE*cos(th3) - l_OB*cos(th3);
                float dxFootL = Jx_th3*dth3 + Jx_th4*dth4;
                float dyFootL = Jy_th3*dth3 + Jy_th4*dth4;

                // Set gains based on buffer and traj times, then calculate desired x,y from Bezier trajectory at current time if necessary
                float teff  = 0;
                float vMult = 0;
                if( t < start_period) { //before sim, waiting on commands
                    if (K_xx > 0 || K_yy > 0) {
                        K_xx = 100; 
                        K_yy = 100; 
                        D_xx = 5;  
                        D_yy = 5;  
                        K_xy = 0;
                        D_xy = 0;
                    }
                    teff = 0;
                }
                else if (t < start_period + traj_period) //within sim time
                {
                    K_xx = input_params[7];  // Foot stiffness N/m
                    K_yy = input_params[8];  // Foot stiffness N/m
                    K_xy = input_params[9];  // Foot stiffness N/m
                    D_xx = input_params[10];  // Foot damping N/(m/s)
                    D_yy = input_params[11];  // Foot damping N/(m/s)
                    D_xy = input_params[12]; // Foot damping N/(m/s)
                    teff = (t-start_period);
                    vMult = 1;
                }
                /*
                else //shuts off
                {
                    teff = traj_period;
                    vMult = 0;
                }
                */
                
                // Get desired foot positions and velocities
                float rDesFoot_L[2] , vDesFoot_L[2];
                float rDesFoot_R[2] , vDesFoot_R[2];
                
                // stuff related to bezier curve foot position/velocity evaluation
                //rDesFoot_bez.evaluate(teff/traj_period,rDesFoot);
                //rDesFoot_bez.evaluateDerivative(teff/traj_period,vDesFoot);
                //vDesFoot[0]/=traj_period;
                //vDesFoot[1]/=traj_period;
                //vDesFoot[0]*=vMult;
                //vDesFoot[1]*=vMult;

                /*
                in order to cycle our leg trajectory don't condition the while to end i.e. use while(t.read())
                then have the while loop be broken based on completion of a certain amount of leg trajectory cycles i.e. 5
                introduce a variable teff=t-start_period and have the desired foot trajectory be based on teff
                use the following for desired foot position/velocity 
                rEd_l = [p_traj.x_0_l; p_traj.y_0_l] + rotation_matrix*p_traj.r_l * [a*cos(omega_swing_l * teff + p_traj.offset); b*sin(omega_swing_l * teff + p_traj.offset)];

                vEd_l = rotation_matrix*p_traj.r_l*[-a*sin(omega_swing_l*teff+ p_traj.offset)*omega_swing_l    ...
                     b*cos(omega_swing_l*teff+ p_traj.offset)*omega_swing_l ]';
                need to get the following parameters from matlab: x_0_l, y_0_l,r_l,a,omega_swing_l,offset,b,theta(for rotation matrix)
                */ 
                
                // if (omega*teff > 6*3.14159) {break;} // allow for 3 cycles then end the experiment

                float c_cos = cos(omega*teff + offset); //where teff = (t-start_period);
                float c_sin  = sin(omega*teff + offset);

                rDesFoot_L[0] = (x_0 + r*(cos(theta)*a*c_cos - sin(theta)*b*c_sin)); //left leg desired pos, vel
                rDesFoot_L[1] = y_0 + r*(sin(theta)*a*c_cos + cos(theta)*b*c_sin);
                vDesFoot_L[0] = (r*omega*(cos(theta)*(-a)*c_sin - sin(theta)*b*c_cos));
                vDesFoot_L[1] = r*omega*(sin(theta)*(-a)*c_sin + cos(theta)*b*c_cos);

                c_cos = cos(omega*teff); //where teff = (t-start_period); 
                c_sin  = sin(omega*teff);

                rDesFoot_R[0] = x_0 + r*(cos(theta)*a*c_cos - sin(theta)*b*c_sin); //same code for right leg, now without offset in c_cos, c_sin
                rDesFoot_R[1] = y_0 + r*(sin(theta)*a*c_cos + cos(theta)*b*c_sin);
                vDesFoot_R[0] = r*omega*(cos(theta)*(-a)*c_sin - sin(theta)*b*c_cos);
                vDesFoot_R[1] = r*omega*(sin(theta)*(-a)*c_sin + cos(theta)*b*c_cos);

                // Calculate the inverse kinematics (joint positions and velocities) for desired joint angles              
                float xFootR_inv = rDesFoot_R[0];
                float yFootR_inv = rDesFoot_R[1];                
                float l_OE_R = sqrt( (pow(xFootR_inv,2) + pow(yFootR_inv,2)) );
                float alpha_R = abs(acos( (pow(l_OE_R,2) - pow(l_AC,2) - pow((l_OB+l_DE),2))/(-2.0f*l_AC*(l_OB+l_DE)) ));
                float th2_des = -(3.14159f - alpha_R); 
                float th1_des = -((3.14159f/2.0f) + atan2(yFootR_inv,xFootR_inv) - abs(asin( (l_AC/l_OE_R)*sin(alpha_R) )));

                float xFootL_inv = rDesFoot_L[0];
                float yFootL_inv = rDesFoot_L[1];                
                float l_OE_L = sqrt( (pow(xFootL_inv,2) + pow(yFootL_inv,2)) );
                float alpha_L = abs(acos( (pow(l_OE_L,2) - pow(l_AC,2) - pow((l_OB+l_DE),2))/(-2.0f*l_AC*(l_OB+l_DE)) ));
                float th4_des = -(3.14159f - alpha_L); 
                float th3_des = -((3.14159f/2.0f) + atan2(yFootL_inv,xFootL_inv) - abs(asin( (l_AC/l_OE_L)*sin(alpha_L) )));
                
                float dd = (Jx_th1*Jy_th2 - Jx_th2*Jy_th1);
                float dth1_des = (1.0f/dd) * (  Jy_th2*vDesFoot_R[0] - Jx_th2*vDesFoot_R[1] );
                float dth2_des = (1.0f/dd) * ( -Jy_th1*vDesFoot_R[0] + Jx_th1*vDesFoot_R[1] );

                float ddL = (Jx_th3*Jy_th4 - Jx_th4*Jy_th3);
                float dth3_des = (1.0f/ddL) * (  Jy_th4*vDesFoot_L[0] - Jx_th4*vDesFoot_L[1] );
                float dth4_des = (1.0f/ddL) * ( -Jy_th3*vDesFoot_L[0] + Jx_th3*vDesFoot_L[1] );
        
                // Calculate error variables
                float e_xR = rDesFoot_R[0] - xFootR;
                float e_yR = rDesFoot_R[1] - yFootR;
                float de_xR = vDesFoot_R[0] - dxFootR;
                float de_yR = vDesFoot_R[1] - dyFootR;

                float e_xL = rDesFoot_L[0] - xFootL;
                float e_yL = rDesFoot_L[1] - yFootL;
                float de_xL = vDesFoot_L[0] - dxFootL;
                float de_yL = vDesFoot_L[1] - dyFootL;
        
                // Calculate virtual force on foot
                // float fx = K_xx*e_x + K_xy*e_y + D_xx*de_x + D_xy*de_y;
                // float fy = K_xy*e_x + K_yy*e_y + D_xy*de_x + D_yy*de_y;
                // float tau1_des = Jx_th1*fx + Jy_th1*fy;
                // float tau2_des = Jx_th2*fx + Jy_th2*fy;

                float fxL = K_xx*e_xL + K_xy*e_yL + D_xx*de_xL + D_xy*de_yL;
                float fyL = K_xy*e_xL + K_yy*e_yL + D_xy*de_xL + D_yy*de_yL;
                float fxR = K_xx*e_xR + K_xy*e_yR + D_xx*de_xR + D_xy*de_yR;
                float fyR = K_xy*e_xR + K_yy*e_yR + D_xy*de_xR + D_yy*de_yR;
                
                float tau1_des = Jx_th1*fxR + Jy_th1*fyR;
                float tau2_des = Jx_th2*fxR + Jy_th2*fyR;
                float tau3_des = Jx_th3*fxL + Jy_th3*fyL;
                float tau4_des = Jx_th4*fxL + Jy_th4*fyL;

                /* control from simulation is currently not working, so it is commented out
                Matrix pVprElL(2,1);
                float pVLx = K_xx*e_xL;
                float pVLy = K_yy*e_yL;
                pVprElL << pVLx     
                        << pVLy ;
                // F_dampl = [-D_x*(vels(1,:) - vEd_l(1)); -D_y*(vels(2,:) - vEd_l(2))];
                Matrix FdampL(2,1);
                float fdLx = -D_xx*de_xL;
                float fdLy = -D_yy*de_yL;
                // Populate matrix
                FdampL.Clear();
                FdampL << fdLx
                        << fdLy ;
                // Populate Jacobian matrix
                JacobianL.Clear();
                JacobianL << Jx_th3 << Jx_th4
                         << Jy_th3 << Jy_th4;
                //taut = -transpose(J_l(1:2,:)) * pVprEl + transpose(J_l(1:2,:))*F_dampl;
                // Calculate tau
                JacobianLT = MatrixMath::Transpose(JacobianL);
                tauL = JacobianLT*pVprElL + JacobianLT*FdampL;
                // Pulling elements from tauL matrix
                float tau3_des = tauL.getNumber(1,1);
                float tau4_des = tauL.getNumber(2,1);

                // same for right foot
                Matrix pVprElR(2,1);
                float pVRx = K_xx*e_xR;
                float pVRy = K_yy*e_yR;
                pVprElR << pVRx     
                        << pVRy ;
                // F_dampl = [-D_x*(vels(1,:) - vEd_l(1)); -D_y*(vels(2,:) - vEd_l(2))];
                Matrix FdampR(2,1);
                float fdRx = -D_xx*de_xR;
                float fdRy = -D_yy*de_yR;
                // Populate matrix
                FdampR.Clear();
                FdampR << fdRx
                        << fdRy ;
                // Populate Jacobian matrix
                JacobianR.Clear();
                JacobianR << Jx_th1 << Jx_th2
                         << Jy_th1 << Jy_th2;
                //taut = -transpose(J_l(1:2,:)) * pVprEl + transpose(J_l(1:2,:))*F_dampl;
                // Calculate tau
                JacobianRT = MatrixMath::Transpose(JacobianR);
                tauR = JacobianRT*pVprElR + JacobianRT*FdampR;
                // Pulling elements from tauL matrix
                float tau1_des = tauR.getNumber(1,1);
                float tau2_des = tauR.getNumber(2,1);
                */

                /*               
                // Calculate mass matrix elements
                float M11 = I1 + I2 + I3 + I4 + Ir + Ir*pow(N,2) + pow(l_AC,2)*m4 + pow(l_A_m3,2)*m3 + pow(l_B_m2,2)*m2 + pow(l_C_m4,2)*m4 
                    + pow(l_OA,2)*m3 + pow(l_OB,2)*m2 + pow(l_OA,2)*m4 + pow(l_O_m1,2)*m1 + 2*l_C_m4*l_OA*m4 + 2*l_AC*l_C_m4*m4*cos(th2) 
                    + 2*l_AC*l_OA*m4*cos(th2) + 2*l_A_m3*l_OA*m3*cos(th2) + 2*l_B_m2*l_OB*m2*cos(th2); 
                float M12 = I2 + I3 + pow(l_AC,2)*m4 + pow(l_A_m3,2)*m3 + pow(l_B_m2,2)*m2 + Ir*N + l_AC*l_C_m4*m4*cos(th2)
                    + l_AC*l_OA*m4*cos(th2) + l_A_m3*l_OA*m3*cos(th2) + l_B_m2*l_OB*m2*cos(th2);
                float M22 = Ir*pow(N,2) + m4*pow(l_AC,2) + m3*pow(l_A_m3,2) + m2*pow(l_B_m2,2) + I2 + I3;
                
                // Populate mass matrix
                MassMatrix.Clear();
                MassMatrix << M11 << M12
                           << M12 << M22;
                
                // Populate Jacobian matrix
                Jacobian.Clear();
                Jacobian << Jx_th1 << Jx_th2
                         << Jy_th1 << Jy_th2;
                
                // Once you have copied the elements of the mass matrix, uncomment the following section
                
                // Calculate Lambda matrix
                JacobianT = MatrixMath::Transpose(Jacobian);
                InverseMassMatrix = MatrixMath::Inv(MassMatrix);
                temp_product = Jacobian*InverseMassMatrix*JacobianT;
                Lambda = MatrixMath::Inv(temp_product); 
                
                // Pull elements of Lambda matrix
                float L11 = Lambda.getNumber(1,1);
                float L12 = Lambda.getNumber(1,2);
                float L21 = Lambda.getNumber(2,1);
                float L22 = Lambda.getNumber(2,2);               
                */

                // float tau1_des = (Jx_th1*L11 + Jx_th2*L21)*fx + (Jx_th1*L12 + Jx_th2*L22)*fy;
                // float tau2_des = (Jy_th1*L11 + Jy_th2*L21)*fx + (Jy_th1*L12 + Jy_th2*L22)*fy;
                                
                // Set desired currents             
                current_des1 = tau1_des/k_t;          
                current_des2 = tau2_des/k_t;
                current_des3 = tau3_des/k_t;
                current_des4 = tau4_des/k_t;

                // Form output to send to MATLAB     
                float output_data[NUM_OUTPUTS];
                // current time
                output_data[0] = t.read();
                // motor 1 state
                output_data[1] = angle1;
                output_data[2] = velocity1;  
                output_data[3] = current1;
                output_data[4] = current_des1;
                output_data[5] = duty_cycle1;
                // motor 2 state
                output_data[6] = angle2;
                output_data[7] = velocity2;
                output_data[8] = current2;
                output_data[9] = current_des2;
                output_data[10]= duty_cycle2;
                //motor 3 state
                output_data[11] = angle3;
                output_data[12] = velocity3;  
                output_data[13] = current3;
                output_data[14] = current_des3;
                output_data[15] = duty_cycle3;
                // motor 4 state
                output_data[16] = angle4;
                output_data[17] = velocity4;
                output_data[18] = current4;
                output_data[19] = current_des4;
                output_data[20]= duty_cycle4;
                // foot state
                output_data[21] = xFootR;
                output_data[22] = yFootR;
                output_data[23] = dxFootR;
                output_data[24] = dyFootR;
                output_data[25] = xFootL;
                output_data[26] = yFootL;
                output_data[27] = dxFootL;
                output_data[28] = dyFootL;
                // desired foot state
                output_data[29] = rDesFoot_R[0];
                output_data[30] = rDesFoot_R[1];
                output_data[31] = vDesFoot_R[0];
                output_data[32] = vDesFoot_R[1];
                output_data[33] = rDesFoot_L[0];
                output_data[34] = rDesFoot_L[1];
                output_data[35] = vDesFoot_L[0];
                output_data[36] = vDesFoot_L[1];
                // desired angles
                output_data[37] = th1_des;
                output_data[38] = dth1_des;

                output_data[39] = th2_des;
                output_data[40] = dth2_des;

                output_data[41] = th3_des;
                output_data[42] = dth3_des;

                output_data[43] = th4_des;
                output_data[44] = dth4_des;
                
                // Send data to MATLAB
                server.sendData(output_data,NUM_OUTPUTS);

                wait_us(impedance_control_period_us);   
            }
            
            // Cleanup after experiment
            server.setExperimentComplete();
            currentLoop.detach();
            motorShield.motorAWrite(0, 0); //turn motor A off
            motorShield.motorBWrite(0, 0); //turn motor B off
            motorShield.motorCWrite(0, 0); //turn motor C off
            motorShield.motorDWrite(0, 0); //turn motor D off
        
        } // end if
        
    } // end while
    
} // end main