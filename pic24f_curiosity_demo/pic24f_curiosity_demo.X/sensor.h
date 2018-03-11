/* 
 * File:   sensor.h
 * Author: Mike
 *
 * Created on January 9, 2018, 6:42 PM
 */

#ifndef SENSOR_H
#define	SENSOR_H
#include "adc.h"
#include <stdint.h>
#include <stdbool.h>
#ifdef	__cplusplus
extern "C" {
#endif

bool SENSOR_Enable();
//sets up the pins and the adc for the sensor

uint16_t SENSOR_ReadValue();
//returns the value from the 

#ifdef	__cplusplus
}
#endif

#endif	/* SENSOR_H */

