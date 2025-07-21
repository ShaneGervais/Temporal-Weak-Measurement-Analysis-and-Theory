from pylablib.devices import Thorlabs #pylablib.devices.Thorlabs.kinesis.BasicKinesisDevice.blink
import pyvisa as visa # http://github.com/hgrecco/pyvisa
import matplotlib.pyplot as plt # http://matplotlib.org/
import numpy as np # http://www.numpy.org/
import time
import csv
from pathlib import Path
import winsound

f = 1000  # Set Frequency
dur = 500  # Set Duration 
duration_meas = [] # I print lists of durations of each measurement and transfer from scope to PC at the end to check 
duration_transfer = []

# Scope sync (OPC)/error register (ESR) function
def scope_sync(type):
    result='Success'
    while True:
        try:
            scope.query('*{}?'.format(type))
            break
        except visa.VisaIOError:
            1#print("Timeout {}".format(type))
    r = scope.query('allev?').strip()
    if not 'No events to report' in r: # For successful measurement the event register either contains 0 or 1 (I don't know the difference)
        winsound.Beep(f, dur) # If it beeps, the measurement was unsuccessful but if automatically takes the same measurement again
        print(r)
        result='Failed'
    scope.write('*cls')  # Clear Error Status Register
    return(result)

# Scope function
def measure_curve():
    scope.write('acquire:state 0')
    scope.write('acquire:stopafter SEQUENCE') # Takes single sequence measurement and then stops
    scope.write('acquire:state 1') # Measurement go
    before_meas = time.time()
    time.sleep(150) # Adjust this to how long a single measurement roughly takes or slightly shorter
    result1 = scope_sync('OPC') # Sync to check that measurement is complete
    after_meas = time.time()
    i=0
    while i<10:
        try:
            waveform =  scope.query_ascii_values('CURVe?', container=np.array, converter='d')
            break
        except visa.VisaIOError:
            print('Problem with data query')
            time.sleep(3)
            scope_sync('OPC')
            scope_sync('ESR')
            print('Sync done')
    result2 = scope_sync('ESR') # Check if data transfer is finished
    after_transfer = time.time()
    duration_meas.append(after_meas-before_meas)
    duration_transfer.append(after_transfer-after_meas)
    #scale waveform
    scaled_time = np.linspace(tstart, tstop, num=np.size(waveform), endpoint=False)
    scaled_wave = (waveform - vpos) * vscale + voff
    if result1 == 'Failed' or result2 == 'Failed':
        return('Failed')
    else:
        return(scaled_time,scaled_wave)

# Rotating mount function
def rotate(stage,target_angle):
    Thorlabs.kinesis.KinesisMotor.move_to(stage,position=target_angle, channel=None, scale='stage')
    Thorlabs.kinesis.KinesisMotor.wait_for_stop(stage)
    if abs(Thorlabs.kinesis.KinesisMotor.get_position(stage,scale='stage')-target_angle) > 0.003:
        print('Further loop needed')
        rotate(stage,target_angle)
    return(Thorlabs.kinesis.KinesisMotor.get_position(stage,scale='stage'))

#Save data function
def save_data_to_csv(filename, xdata, ydata, meas_deg):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Measured deg = {}'.format(meas_deg)])
        writer.writerow(['Times/s', 'Voltage/V'])  # Write header
        for j in range(len(ydata)):
            writer.writerow([xdata[j], ydata[j]])


## Setting up stages
# Stages are numbered as they appear in setup, laser passes stage1 first..
stage1=Thorlabs.kinesis.KinesisMotor("27267156", scale="stage", is_rack_system=False)
stage2=Thorlabs.kinesis.KinesisMotor("27267210", scale="stage", is_rack_system=False)
stage3=Thorlabs.kinesis.KinesisMotor("27267041", scale="stage", is_rack_system=False)
stages= [stage1,stage2,stage3]
for stage in stages:
    assert Thorlabs.kinesis.KinesisMotor.get_scale_units(stage) =='deg' # Assert input is understood as physical degrees
    if Thorlabs.kinesis.KinesisMotor.is_homed(stage) == False: # Home stages if needed
        print('stage {} homed successfully.'.format(stage))
        Thorlabs.kinesis.KinesisMotor.home(stage)
        Thorlabs.kinesis.KinesisMotor.wait_for_home(stage)

## Define global variables
Nr_averages_scope = 10000
Nr_averages_pc = 5
rotate(stage2,355) # Manually rotate fixed stages 
rotate(stage3,330)
variable_stage = stage1
variable_stage_offset = 3
min_deg = 0
max_deg = 90
step_nr_deg = 37 # How many steps

## Setting up scope
input("Check that Input is on CH1 and trigger on AUX. Press Enter to continue...")
# Better keep the following
visa_address = 'GPIB0::1::INSTR'  # Find this in NI MAX application, oscilloscope programm must be active
rm = visa.ResourceManager()
scope = rm.open_resource(visa_address) # If you get an error here its bad, remove hardware conection to PC and connect again
scope.timeout = 10000  # ms
scope.write('*cls')  # Clear Error Status Register
scope.write('DATa:ENCdg ASCII') # ASCII workes for now, I had problems with binary
scope.write('wfmoutpre:byt_nr 8')
scope.write('*rst')  # Reset

# Scope Settings that can be changed
scope.write('DATa:SOUrce CH1')
scope.write('CH1:TERmination 50.0E+0') # impedance
scope.write('ACQuire:SAMPlingmode ET') # equivalent time
scope.write('ACQuire:MODe AVErage')
scope.write('acquire:stopafter RUNSTop') # single
scope.write('ACQuire:NUMAVg {}'.format(Nr_averages_scope))  # Set the number of waveform acquisitions to average
scope.write('CH1:scale 5e-3') # in V
scope.write('CH1:position -3') # in sections (not V)
scope.write('HORizontal:RESOlution 5000') # this is the 'Rec Length' number on the scope
scope.write('horizontal:scale 2e-9')
scope.write('DATa:STARt 1') # data point nr 1 is the start
record_length = int(scope.query('horizontal:recordlength?')) 
scope.write('DATa:STOP {}'.format(record_length)) # last data point is end
# Trigger
scope.write('TRIGger:A:EDGE:SOUrce AUXiliary')
scope.write('TRIGger:A:EDGE:COUPling DC')
scope.write('TRIGger:A:LEVel 0.370')
scope.write('HORizontal:DELay:MODe ON')
scope.write('HORizontal:DELay:POSition 20')
print('Scope settings updated successfully.')

# Retrieve scope settings for scaling
tscale = float(scope.query('wfmoutpre:xincr?'))
tstart = float(scope.query('wfmoutpre:xzero?'))
vscale = float(scope.query('wfmoutpre:ymult?'))  # volts / level
voff = float(scope.query('wfmoutpre:yzero?'))  # reference voltage
vpos = float(scope.query('wfmoutpre:yoff?'))  # reference position (level)
total_time = tscale * record_length
tstop = tstart + total_time    

# Create folder for measurement
folder_name = input('Enter folder name: ')
if folder_name == '':
    print("Folder name cannot be empty. Try again")
    folder_name = input('Enter folder name: ')

folder_path = Path(r'C:\Users\aew8076\Documents\Good measurements') / folder_name

# If folder exists, clear existing contents. Otherwise, create it
if folder_path.exists():
    input("Overwriting old folder. Press Enter to continue...")
    for item in folder_path.iterdir():
        item.unlink()
else:
    folder_path.mkdir(parents=True, exist_ok=True)

# Perform measurement
scope_sync('OPC')
degrees = np.linspace(min_deg+variable_stage_offset,max_deg+variable_stage_offset,step_nr_deg)
print(degrees)
times_list = []
scope.write('acquire:state 0') # Stop scope
for deg in degrees:
    print('Starting {}°'.format(deg))
    averaged_wave = np.zeros(record_length)
    meas_deg=rotate(variable_stage,deg)
    for i in range(Nr_averages_pc):
        while True:
            start=time.time()
            # print('Starting measurement Nr.', i+1)
            measurements = measure_curve()
            if measurements == 'Failed':
                print('Measurement for {} failed! Retrying measurement {}...'.format(deg,i))
            else:
                times,wave = measurements
                break  # Break out of the retry loop if measurement succeeded
        # Save data to CSV file
        filename = folder_path / f'mesure_faible_HeNe_5cm_{int(deg)}_deg_{i+1}.csv'
        save_data_to_csv(filename, times, wave, meas_deg)
        averaged_wave += wave
        end=time.time()
        duration= end-start
        times_list.append(duration)
        scope_sync('OPC')
        if deg==min_deg+variable_stage_offset and i==0:
            print('Estimated time in min: ', (duration*Nr_averages_pc*step_nr_deg + step_nr_deg*2 - duration)/60)
    averaged_wave /= Nr_averages_pc
    #save_data_to_csv(folder_path / f'mesure_faible_HeNe_5cm_{deg}_deg_avg.csv',times,averaged_wave,meas_deg)
    print('Finished {}°'.format(deg))

print('Measurement complete, go to {}'.format(folder_name))

print('Duration of single measurements',duration_meas)
print('Duration of data transfer',duration_transfer)

# Close stages
stage1.close()
stage2.close()
stage3.close()

# Check error status and close scope
while True:
    try:
        r = int(scope.query('*esr?'))
        break
    except visa.VisaIOError:
        print('Timeout while checking error status 1')
print('Event status register: 0b{:08b}'.format(r))
while True:
    try:
        r = scope.query('allev?').strip()
        break
    except visa.VisaIOError:
        print('Timeout while checking error status 2')
print('All event messages: {}'.format(r))
scope.close()
rm.close()