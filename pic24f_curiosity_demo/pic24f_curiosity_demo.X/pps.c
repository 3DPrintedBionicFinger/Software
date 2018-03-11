#include <xc.h>
#include "pps.h"
#include <stdint.h>
#include <stdbool.h>

#define PIN_ANALOG      1
#define PIN_DIGITAL     0

#define PIN_INPUT       1
#define PIN_OUTPUT      0

 bool PPS_OutputSelect(OUTPUT_SOURCE device, PIN_MAPPED pin, bool override){
     switch(pin)
     {
        case RP0:
            if(!(RPOR0bits.RP0R == 0)){
                if(!override)
                    return true;
                RPOR0bits.RP0R=device;
                return true;
            }
            RPOR0bits.RP0R=device;
            break;
        case RP1:
            if(!(RPOR0bits.RP1R == 0)){
                if(!override)
                    return true;
                RPOR0bits.RP1R=device;
                return true;
            }
            RPOR0bits.RP1R=device;
            break;
        case RP2:
            if(!(RPOR1bits.RP2R == 0)){
                if(!override)
                    return true;
                RPOR1bits.RP2R=device;
                return true;
            }
            RPOR1bits.RP2R=device;
            break;
        case RP3:
            if(!(RPOR1bits.RP3R == 0)){
                if(!override)
                    return true;
                RPOR1bits.RP3R=device;
                return true;
            }
            RPOR1bits.RP3R=device;
            break;
        case RP5:
            if(!(RPOR2bits.RP5R == 0)){
                if(!override)
                    return true;
                RPOR2bits.RP5R=device;
                return true;
            }
            RPOR2bits.RP5R=device;
            break;
        case RP6:
            if(!(RPOR3bits.RP6R == 0)){
                if(!override)
                    return true;
                RPOR3bits.RP6R=device;
                return true;
            }
            RPOR3bits.RP6R=device;
            break;
        case RP7:
            if(!(RPOR3bits.RP7R == 0)){
                if(!override)
                    return true;
                RPOR3bits.RP7R=device;
                return true;
            }
            RPOR3bits.RP7R=device;
            break;
        case RP8:
            if(!(RPOR4bits.RP8R == 0)){
                if(!override)
                    return true;
                RPOR4bits.RP8R=device;
                return true;
            }
            RPOR4bits.RP8R=device;
            break;
        case RP9:
            if(!(RPOR4bits.RP9R == 0)){
                if(!override)
                    return true;
                RPOR4bits.RP9R=device;
                return true;
            }
            RPOR4bits.RP9R=device;
            break;
        case RP10:
            if(!(RPOR5bits.RP10R == 0)){
                if(!override)
                    return true;
                RPOR5bits.RP10R=device;
                return true;
            }
            RPOR5bits.RP10R=device;
            break;
        case RP11:
            if(!(RPOR5bits.RP11R == 0)){
                if(!override)
                    return true;
                RPOR5bits.RP11R=device;
                return true;
            }
            RPOR5bits.RP11R=device;
            break;
     }
     return true;
 }