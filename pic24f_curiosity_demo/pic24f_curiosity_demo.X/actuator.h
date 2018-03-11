/* 
 * File:   actuator.h
 * Author: Mike
 *
 * Created on January 12, 2018, 11:54 AM
 */

#ifndef ACTUATOR_H
#define	ACTUATOR_H

#include <stdint.h>
#include <stdbool.h>

#ifdef	__cplusplus
extern "C" {
#endif

bool ACTUATOR_Enable();
//sets up the pins and the adc for the sensor

bool ACTUATOR_WriteValue();
//returns the value from the 


#ifdef	__cplusplus
}
#endif

#endif	/* ACTUATOR_H */

