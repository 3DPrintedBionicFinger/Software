#include <xc.h>
#include "sensor.h"
#include <stdint.h>
#include <stdbool.h>

#define PIN_ANALOG      1
#define PIN_DIGITAL     0

#define PIN_INPUT       1
#define PIN_OUTPUT      0

#define ADC_CHANNEL_SENSOR ADC_CHANNEL_8 //RB12
#define SENSOR_TRIS TRISBbits.TRISB12

bool SENSOR_Enable(){
    SENSOR_TRIS = PIN_INPUT;
    return (ADC_ChannelEnable(ADC_CHANNEL_SENSOR));
}

uint16_t SENSOR_ReadValue(){
    uint16_t temp = 0;
    temp = ADC_Read10bit(ADC_CHANNEL_SENSOR);
    
    //Linerize data here
    
    return temp;
}