#define SEMG_PIN A0
#define SENSOR_PIN A1
#define ACTUATOR_PIN 11
#define H_BRIDGE_FORWARD_PIN
#define H_BRIDGE_BACKWARD_PIN

//may set up a calibration at the start
#define REST_FREQUECY_BUTTON_PIN 0
#define MAX_FORCE_FREQUECY_BUTTON_PIN 1
#define CALIBRATED_BUTTON_PIN 2
#define BUTTON_PUSHED LOW


//all serail prints happen within the ifndef statments
#define DEBUG//needed for the others, dont have more then one on at the same time, serial will be a mess of values
//#define DEBUG_ACT
#define DEBUG_SEMG
//#define DEBUG_SENSOR
//#define DEBUG_TASK_TIMES//messes with all timing and interupts when on
//#define DEBUG_FFT


//#define CALIBRATE //not sure how it reacts when the circuit is not set up



int sliceNum=0;

#define SAMPLE_RATE 500//
#define BUFFER_SIZE 128//make sure to change BUFFER_SIZE_POWER
#define BUFFER_SIZE_POWER 7 //make sure to change BUFFER_SIZE
int FFTBuffer[BUFFER_SIZE]={0};
int bufferCount=0;

float FFTx[BUFFER_SIZE]={0};
float FFTy[BUFFER_SIZE]={0};
float frequency=60;

int sensorBuffer=0;
#define SENSOR_MIN 0//define limits of the feed back sensor that could break the actuator, used for testing, thats all
#define SENSOR_MAX 100

#define POSITION_MAX SENSOR_MIN //values of the feedback sensor that could break the finger, will be different from SENSOR_MIN and MAX
#define POSITION_MIN SENSOR_MAX
#define POSITION_BUFFER 5//adds a value to the min, and subtracts from the max


#define FORCE_CONTROL


#define DUTY_CYCLE_MAX 255 //the max force we want to apply to prevent breaking
#define DUTY_CYCLE_REST 100 // the 
#define DUTY_CYCLE_MIN 0

int actDutyCycle=0;

bool calibrated=false;
float restFrequency=60;
float maxForceFrequency=160;

float kP=1, kI=0, kD=0;
float error1=0;
float error2=0;

#ifdef DEBUG_TASK_TIMES
unsigned long readSEMGTime=0, writeActTime=0, FFTTime=0, controlTime=0;
#endif

void readSEMG(void);
void readSensor(void);
void writeActuator(void);



void setup() {
  #ifdef CALIBRATE
  pinMode(REST_FREQUECY_BUTTON_PIN, INPUT);
  pinMode(MAX_FORCE_FREQUECY_BUTTON_PIN, INPUT);
  pinMode(CALIBRATED_BUTTON_PIN, INPUT);
  #endif
  // put your setup code here, to run once:
  cli();//stop interrupts
  //set timer1 interrupt at 1KHz
  TCCR1A = 0;// set entire TCCR1A register to 0
  TCCR1B = 0;// same for TCCR1B
  TCNT1  = 0;//initialize counter value to 0
  // set compare match register for 1hz increments
  OCR1A = 15999;// = (16*10^6) - 1 (must be <65536) 15999
  // turn on CTC mode
  TCCR1B |= (1 << WGM12);
  // Set CS10 bit for no prescaler
  TCCR1B |= (1 << CS10);  
  // enable timer compare interrupt
  TIMSK1 |= (1 << OCIE1A);
  
  #ifdef DEBUG
  Serial.begin(9600);
  Serial.print("SEMG_PIN=");
  Serial.print(SEMG_PIN);
  Serial.print("\nSENSOR_PIN=");
  Serial.print(SENSOR_PIN);
  Serial.print("\nACTUATOR_PIN=");
  Serial.print(ACTUATOR_PIN);
  Serial.print("\n");
  #endif
  #ifdef DEBUG_ACT
  Serial.print("actuatorDutyCycle\n");
  #endif
  #ifdef DEBUG_SEMG
  Serial.print("FFTBuffer,bufferCount\n");
  #endif
  #ifdef DEBUG_SENSOR
  Serial.print("sensorBuffer\n");
  #endif
  #ifdef DEBUG_FFT
  Serial.print("FFT frequency\n");
  #endif
  #ifdef DEBUG_TASK_TIMES
  Serial.print("Read SEMG Time,  Actuator Write Time, FFT Time,  Control Time\n");
  #endif
//  #ifdef DEBUG_FFT
  
  sei();//allow interrupts
  delay(250);//this is to allow the FFTBuffer to fill for the first time  
}

void loop() {
  // put your main code here, to run repeatedly:
  #ifdef DEBUG_TASK_TIMES
  cli();
  FFTTime=micros();
  #endif
  int FFTEnd=bufferCount;//this is to prevent inturupps from messing with the order while filling FFTx
  for(int i =0; i<BUFFER_SIZE; i++){   //rearanges the buffer to proper order
    int index=(i+FFTEnd+1)%BUFFER_SIZE;
    FFTx[i]=FFTBuffer[index];
    FFTy[i]=0;
  }
  FFT(1,BUFFER_SIZE_POWER,FFTx,FFTy);
  frequency=FFTFrequency();
  #ifdef DEBUG_FFT
  Serial.print(frequency);
  Serial.print("\n");
  #endif
  //frequency=goertzel_mag(BUFFER_SIZE,60.0,SAMPLE_RATE, FFTx);//float goertzel_mag(int numSamples,float TARGET_FREQUENCY,int SAMPLING_RATE, float* data)
  #ifdef DEBUG_TASK_TIMES
  FFTTime=micros()-FFTTime;
  Serial.print(readSEMGTime);
  Serial.print(",");
  Serial.print(writeActTime);
  Serial.print(",");
  Serial.print(FFTTime);
  Serial.print(",");
  Serial.print(controlTime);
  Serial.print("\n");
  sei();
  delay(250);
  #endif
}

ISR(TIMER1_COMPA_vect){//timer1 interrupt 1KHz 
  switch(sliceNum){//this is a time slice schedular, with 1ms slices
    case 0:
      #ifdef DEBUG_TASK_TIMES
      cli();
      readSEMGTime=micros();
      #endif
      readSEMG();
      #ifdef DEBUG_TASK_TIMES
      readSEMGTime=micros()-readSEMGTime;
      sei();
      #endif
      
      break;
     case 1:
      readSensor();
      break;
     case 2:
      readSEMG();
      break;
     case 3:
     #ifdef DEBUG_TASK_TIMES
      cli();
      controlTime=micros();
      #endif
      runControl();
      #ifdef DEBUG_TASK_TIMES
      controlTime=micros()-controlTime;
      sei();
      #endif
      break;
     case 4:
      readSEMG();
      break;
     case 5:
     #ifdef DEBUG_TASK_TIMES
      cli();
      writeActTime=micros();
      #endif
      writeActuator();
      #ifdef DEBUG_TASK_TIMES
      writeActTime=micros()-writeActTime;
      sei();
      #endif
      break;
     case 6:
      readSEMG();
      break;
     case 7:
      #ifdef CALIBRATE
      if(digitalRead(CALIBRATED_BUTTON_PIN)==BUTTON_PUSHED){
          if(calibrated){
            calibrated=false;
          }else{
            calibrated=true;
          }
        }
      if(calibrated!){
        if(digitalRead(REST_FREQUECY_BUTTON_PIN)==BUTTON_PUSHED){
          restFrequency=frequency;
        }
        if(digitalRead(MAX_FORCE_FREQUECY_BUTTON_PIN)==BUTTON_PUSHED){
          maxForceFrequency=frequency;
        }
      }
      #endif
     case 8:
      readSEMG();
      break;
     case 10:
      readSEMG();
      break;
     case 12:
      readSEMG();
      break;
     case 14:
      readSEMG();
      break;
     case 16:
      readSEMG();
      break;
     case 18:
      readSEMG();
      break;
     case 20:
      readSEMG();
      break;
     case 22:
      readSEMG();
      break;
     case 24:
      readSEMG();
      break;
     case 26:
      readSEMG();
      break;
     case 28:
      readSEMG();
      break;
     default:
      break;
  }
  sliceNum++;
  sliceNum=sliceNum%30;
}
void readSEMG(void){
  FFTBuffer[bufferCount]=analogRead(SEMG_PIN);
  #ifdef DEBUG_SEMG
  Serial.print(FFTBuffer[bufferCount]);
  Serial.print(",");
  Serial.print(bufferCount);
  Serial.print("\n");
  #endif
  bufferCount++;
  bufferCount=bufferCount%BUFFER_SIZE;
  
}
void readSensor(void){
  sensorBuffer=analogRead(SENSOR_PIN);
  #ifdef DEBUG_SENSOR
  Serial.print(sensorBuffer);
  //Serial.print(",");
  Serial.print("\n");
  #endif
}
void writeActuator(void){
  analogWrite(ACTUATOR_PIN,actDutyCycle);
  #ifdef DEBUG_ACT
  Serial.print(actDutyCycle);
  //Serial.print(",");
  Serial.print("\n");
  #endif
}
void runControl(){
  float output=0;
  float input=0;
  float sensor=0;
  input=frequencyLinerize(frequency);
  sensor=sensorLinerize(sensorBuffer);
  error2=error1;
  error1=input-sensor;
  output=kP*error1+kI*(error1+error2)+kD*(error1-error2);
  
  #ifdef FORCE_CONTROL
  
  #endif
}

short FFT(short int dir,int m,float *x,float *y)
{
   long n,i,i1,j,k,i2,l,l1,l2;
   float c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
   
   return(true);
}
float FFTFrequency(){
  int maxIndex=0;
  float maxMagnitude=0;
  float currentMagnitude=0;
  for(int i=4; i<(BUFFER_SIZE/2);i++){//start at 4 is to skip problems at the lower frequncies due to descretization
    currentMagnitude=FFTx[i]*FFTx[i]+FFTy[i]*FFTy[i];
    if(currentMagnitude>maxMagnitude){
      maxIndex=i;
      maxMagnitude=currentMagnitude;
    }
  }
  return maxIndex*SAMPLE_RATE/(BUFFER_SIZE/2);
  
}
float frequencyLinerize(float frequency){//gives a value from 0-1 of frequency
  float output=0;
  if (frequency < restFrequency){
    frequency=restFrequency;
  }
  if (frequency > maxForceFrequency){
    frequency=maxForceFrequency;
  }

  output=(frequency-restFrequency)/(maxForceFrequency-restFrequency);
  return output;
  
  
}

float sensorLinerize(int input){//gives a value from 0-1 based on senor value
  float output=0; 
  
  output=(input-(POSITION_MIN+POSITION_BUFFER))/(POSITION_MAX-POSITION_MIN-POSITION_BUFFER*2);
  return output;
}

int controlToActDutyCycle(float input){
  
}
float goertzel_mag(int numSamples,float TARGET_FREQUENCY,int SAMPLING_RATE, float* data)
{
    int     k,i;
    float   floatnumSamples;
    float   omega,sine,cosine,coeff,q0,q1,q2,magnitude,real,imag;

    float   scalingFactor = numSamples / 2.0;

    floatnumSamples = (float) numSamples;
    k = (int) (0.5 + ((floatnumSamples * TARGET_FREQUENCY) / (float)SAMPLING_RATE));
    omega = (2.0 * M_PI * k) / floatnumSamples;
    sine = sin(omega);
    cosine = cos(omega);
    coeff = 2.0 * cosine;
    q0=0;
    q1=0;
    q2=0;

    for(i=0; i<numSamples; i++)
    {
        q0 = coeff * q1 - q2 + data[i];
        q2 = q1;
        q1 = q0;
    }

    // calculate the real and imaginary results
    // scaling appropriately
    real = (q1 - q2 * cosine) / scalingFactor;
    imag = (q2 * sine) / scalingFactor;

    magnitude = sqrtf(real*real + imag*imag);
    return magnitude;
}

