/*******************************************************************************
Copyright 2016 Microchip Technology Inc. (www.microchip.com)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*******************************************************************************/

#include <stdint.h>

#include "leds.h"
#include "buttons.h"
#include "leds_rgb.h"
#include "adc.h"
#include "timer_1ms.h"
#include "sensor.h"

//------------------------------------------------------------------------------
//Application related definitions
//------------------------------------------------------------------------------
#define BUTTON_DEBOUCE_TIME_MS      20

typedef enum
{
    BUTTON_COLOR_RED = 0,
    BUTTON_COLOR_GREEN = 1,
    BUTTON_COLOR_BLUE = 2
} BUTTON_COLOR;

#define MS_PER_TICK 100
#define TICK_MAX 60
#define TICKS_PER_SMEG_INPUT 1
#define TICKS_PER_SENSOR_INPUT 5
#define TICKS_PER_ACTUATOR_OUTPUT TICKS_PER_SENSOR_INPUT
#define TICKS_PER_CONTROL_CALC TICKS_PER_SENSOR_INPUT
#define TICKS_PER_FFT TICKS_PER_SENSOR_INPUT

//------------------------------------------------------------------------------
//Private prototypes
//------------------------------------------------------------------------------
static void ButtonDebounce(void);
static void TimerISR(void);

//------------------------------------------------------------------------------
//Global variables
//------------------------------------------------------------------------------
static volatile BUTTON_COLOR buttonColor;
static volatile int tickCount = 0;
static volatile int mSCount = 0;

//these are for running tasks
static volatile bool feedBackTask = false;
static volatile bool outputTask = false;
static volatile bool runControlTask = false;
static volatile bool fFTTask = false;


//------------------------------------------------------------------------------
//Main Function
//------------------------------------------------------------------------------
int main(void)
{
    /*
    uint16_t potentiometer;
    uint16_t red = 512, green = 256, blue = 100;
    
    buttonColor = BUTTON_COLOR_RED;
    
    //Configure the pushbutton pins as digital inputs.
    BUTTON_Enable(BUTTON_S1);
    BUTTON_Enable(BUTTON_S2);

    //Configure and enable the I/O pins controlling the general purpose LEDs, and 
    //the PWM outputs controlling the RGB LED.
    LED_Enable(LED_LD1);
    LED_Enable(LED_LD2);
    LED_RGB_Enable(LED_RGB_LD3);
    
    //Enable and configure the ADC so it can sample the potentiometer.
    ADC_SetConfiguration(ADC_CONFIGURATION_DEFAULT);
    ADC_ChannelEnable(ADC_CHANNEL_POTENTIOMETER);
    
    //Turn on a timer, so to generate periodic interrupts.
    TIMER_SetConfiguration(TIMER_CONFIGURATION_1MS);
    
    //Register the ButtonDebounce() callback function, so it gets called periodically
    //when the timer interrupts occur (in this case at 1:1 rate, so ButtonDebounce()
    //executes once per 1ms).
    TIMER_RequestTick(&ButtonDebounce, 1);
    */
    /***********************************************************************/
    
    
    //Set the default config on the adc.
    ADC_SetConfiguration(ADC_CONFIGURATION_DEFAULT);
    
    //set up the sensor, needs adc set up first(probably).
    SENSOR_Enable();
    
    //Turn on a timer, so to generate periodic interrupts.
    TIMER_SetConfiguration(TIMER_CONFIGURATION_1MS);
    
    //I think this is used to register the function called by the ISR for the timer
    TIMER_RequestTick(&TimerISR, 1);
    //Main loop
    while(1)
    {
        if (feedBackTask){
            feedBackTask=false;
            
        }else if (runControlTask){
            runControlTask = false;
            
        }else if (outputTask){
            outputTask=false;
            
        }else if (fFTTask){
            fFTTask=false;
            
        }
        
        /**********************************************************************/
        /*
        //Fetch an ADC sample from the potentiometer
        potentiometer = ADC_Read10bit(ADC_CHANNEL_POTENTIOMETER);

        //Use the potentiometer ADC value to set the brightness of the currently
        //selected color channel on the RGB LED.  The "currently selected channel"
        //is manually selected by the user at runtime by pressing the pushbuttons.
        switch(buttonColor)
        {
            case BUTTON_COLOR_RED:
                red = potentiometer;
                break;

            case BUTTON_COLOR_GREEN:
                green = potentiometer;
                break;

            case BUTTON_COLOR_BLUE:
                blue = potentiometer;
                break;

            default:
                break;
        }

        //Update the PWM values controlling the intensity of the RGB LED channels.
        LED_RGB_Set(LED_RGB_LD3, red, green, blue);
        */
        
    }
}


//Helper function that advances the currently selected RGB color channel that
//is to be adjusted next.  This function is called in response to user pushbutton
//press events.
static void ChangeColor(void)
{         
    switch(buttonColor)
    {
        case BUTTON_COLOR_RED:
            buttonColor = BUTTON_COLOR_GREEN;
            break;

        case BUTTON_COLOR_GREEN:
            buttonColor = BUTTON_COLOR_BLUE;
            break;

        case BUTTON_COLOR_BLUE:
            buttonColor = BUTTON_COLOR_RED;
            break;

        default:
            buttonColor = BUTTON_COLOR_RED;
            break;
    }
}


//This callback function gets called periodically (1/1ms) by the timer interrupt event
//handler.  This function is used to periodically sample the pushbuttons and implement
//a debounce algorithm to reject spurious chatter that can occur during press events.
static void ButtonDebounce(void)
{
    static uint16_t debounceCounterS1 = 0;
    static uint16_t debounceCounterS2 = 0;
    
    //Sample the button S1 to see if it is currently pressed or not.
    if(BUTTON_IsPressed(BUTTON_S1))
    {
        //The button is currently pressed.  Turn on the general purpose LED.
        LED_On(LED_LD1);
        
        //Check if the debounce blanking interval has been satisfied.  If so,
        //advance the RGB color channel user control selector.
        if(debounceCounterS1 == 0)
        {
            ChangeColor();   
        }
        
        //Reset the debounce countdown timer, so a new color change operation
        //won't occur until the button is released and remains continuously released 
        //for at least BUTTON_DEBOUCE_TIME_MS.
        debounceCounterS1 = BUTTON_DEBOUCE_TIME_MS;
    }
    else
    {
        //The button is not currently pressed.  Turn off the LED.
        LED_Off(LED_LD1);  
        
        //Allow the debounce interval timer to count down, until it reaches 0.
        //Once it reaches 0, the button is effectively "re-armed".
        if(debounceCounterS1 != 0)
        {
            debounceCounterS1--;
        }
    }
    
    
    
    //Sample the button S2 to see if it is currently pressed or not.
    if(BUTTON_IsPressed(BUTTON_S2))
    {
        //The button is currently pressed.  Turn on the general purpose LED.
        LED_On(LED_LD2);
        
        //Check if the debounce blanking interval has been satisfied.  If so,
        //advance the RGB color channel user control selector.
        if(debounceCounterS2 == 0)
        {
            ChangeColor();   
        }
        
        //Reset the debounce countdown timer, so a new color change operation
        //won't occur until the button is released and remains continuously released 
        //for at least BUTTON_DEBOUCE_TIME_MS.
        debounceCounterS2 = BUTTON_DEBOUCE_TIME_MS;
    }
    else
    {
        //The button is not currently pressed.  Turn off the LED.
        LED_Off(LED_LD2); 
        
        //Allow the debounce interval timer to count down, until it reaches 0.
        //Once it reaches 0, the button is effectively "re-armed".
        if(debounceCounterS2 != 0)
        {
            debounceCounterS2--;
        }
    }    
    
}

static void TimerISR(void){
    mSCount = mSCount + 1;
    if (mSCount > MS_PER_TICK ) {
        tickCount++
        mSCount = 0;
        
    }
    
    tickCount = tickCount%TICK_MAX;
    //if tickCount is an even multiple of TICKS_PER_SMEG_INPUT
    if(!(tickCount % TICKS_PER_SMEG_INPUT)){
        //given preemptive abbiltiy for this task by including in the ISR
        
    }
     //if tickCount is an even multiple of TICKS_PER_ACTUATOR_OUTPUT
    if(!(tickCount % TICKS_PER_ACTUATOR_OUTPUT)){
        outputTask=true;
    }
     //if tickCount is an even multiple of TICKS_PER_CONTROL_CALC
    if(!(tickCount % TICKS_PER_CONTROL_CALC)){
        runControlTask=true;
    }
     //if tickCount is an even multiple of TICKS_PER_SENSOR_INPUT
    if(!(tickCount % TICKS_PER_SENSOR_INPUT)){
        feedBackTask=true;
    }
    //if tickCount is an even multiple of TICKS_PER_FFT
    if(!(tickCount % TICKS_PER_FFT)){
        fFTTask=true;
    }
    return;
}