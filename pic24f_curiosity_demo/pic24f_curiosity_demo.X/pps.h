/* 
 * File:   pps.h
 * Author: Mike
 *
 * Created on January 12, 2018, 12:31 PM
 */

#ifndef PPS_H
#define	PPS_H

#include <stdint.h>
#include <stdbool.h>

#ifdef	__cplusplus
extern "C" {
#endif
    
    typedef enum
    {
        OC1=13,
        OC2=14
        
    } OUTPUT_SOURCE;
    
    typedef enum
    {
        RP0,//RB0-Pin 21
        RP1,//RB1-Pin 22
        RP2,//RB2-Pin 23
        RP3,//RB3-Pin 24
        //RP4 is not implemented
        RP5,//RB5-Pin 41 
        RP6,//RB6-Pin 42
        RP7,//RB7-Pin 43
        RP8,//RB8-Pin 44
        RP9,//RB9
        RP10,//RB10
        RP11,//RB11
        //RP12 is not implemented
        RP13,//RB13
        RP14,//RB14
        RP15,//RB15
        RP16//RC0
    } PIN_MAPPED;
    
    /*
     * bool PPS_OutputSelect(OUTPUT_SOURCE, PIN_MAPPED,bool)
     * 
     * OUTPUT_SOURCE device, variable controlling the output device to be connected
     * PIN_MAPPED pin, variable controlling the pin the output device connects to
     * bool override, variable controlling if the device should overide a device already being used
     * 
     * outputs true if the pin selected is already assigned a device,false otherwise
     * 
     * 
     * if override is set to false then the device will not be mapped to the pin
     * if the pin is already in use
     * 
     * if override is set to true, then the device will be mapped to the pin
     * regardless of if the pin is in use
     */
    bool PPS_OutputSelect(OUTPUT_SOURCE device, PIN_MAPPED pin, bool override);
    
    

#ifdef	__cplusplus
}
#endif

#endif	/* PPS_H */

