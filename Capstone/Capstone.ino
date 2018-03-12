#define SEMG_PIN A0
#define SENSOR_PIN A1
#define ACTUATOR_PIN 11



#define DEBUG




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
#define SENSOR_MIN
#define SENSOR_MAX
#define 

int actuatorDutyCycle=0;

#define MAX_POSITION
#define MIN_POSITION

#define FORCE_CONTROL

bool calibrated=false;
float restFrequency=60;
float maxFlexFrequency=200;

float kP=0, kI=0, kD=0;
float error1=0;
float error2=0;

void readSEMG(void);
void readSensor(void);
void writeActuator(void);



void setup() {
  // put your setup code here, to run once:
  cli();//stop interrupts
  //set timer1 interrupt at 1KHz
  TCCR1A = 0;// set entire TCCR1A register to 0
  TCCR1B = 0;// same for TCCR1B
  TCNT1  = 0;//initialize counter value to 0
  // set compare match register for 1hz increments
  OCR1A = 15999;// = (16*10^6) - 1 (must be <65536)
  // turn on CTC mode
  TCCR1B |= (1 << WGM12);
  // Set CS10 bit for no prescaler
  TCCR1B |= (1 << CS12) | (1 << CS10);  
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
  #endif
  sei();//allow interrupts
  delay(250);//this is to allow the FFTBuffer to fill for the first time
}

void loop() {
  // put your main code here, to run repeatedly:
  int FFTEnd=bufferCount;//this is to prevent inturupps from messing with the order while filling FFTx
  for(int i =0; i<BUFFER_SIZE; i++){   //rearanges the buffer to proper order
    int index=(i+FFTEnd+1)%BUFFER_SIZE;
    FFTx[i]=FFTBuffer[index];
    FFTy[i]=0;
  }
  FFT(1,BUFFER_SIZE_POWER,FFTx,FFTy);
  frequency=FFTFrequency();
}

ISR(TIMER1_COMPA_vect){//timer1 interrupt 1KHz 
  switch(sliceNum){
    case 0:
      readSEMG();
      break;
     case 1:
      readSensor();
      break;
     case 2:
      readSEMG();
      break;
     case 3:
      //runControl();
      break;
     case 4:
      readSEMG();
      break;
     case 5:
      writeActuator();
      break;
     case 6:
      readSEMG();
      break;
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
  bufferCount++;
  bufferCount=bufferCount%512;
}
void readSensor(void){
  sensorBuffer=analogRead(SENSOR_PIN);
}
void writeActuator(void){
  analogWrite(ACTUATOR_PIN,actuatorDutyCycle);
}
void runControl(){
  float output=0;

  
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
  for(int i=4; i<(BUFFER_SIZE/2);i++){
    currentMagnitude=FFTx[i]*FFTx[i]+FFTy[i]*FFTy[i];
    if(currentMagnitude>maxMagnitude){
      maxIndex=i;
      maxMagnitude=currentMagnitude;
    }
  }
  return maxIndex*SAMPLE_RATE/(BUFFER_SIZE/2);
  
}
