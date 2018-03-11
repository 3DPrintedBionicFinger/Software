/* 
 * File:   pwm.h
 * Author: Mike
 *
 * Created on January 12, 2018, 11:57 AM
 */

#ifndef PWM_H
#define	PWM_H

#ifdef	__cplusplus
extern "C" {
#endif
    
    typedef enum
    {
        OC1=
    } OCX;

    bool PWM_Configure();
    
    bool PWM_PinSelect(OCX, RPN);

#ifdef	__cplusplus
}
#endif

#endif	/* PWM_H */

