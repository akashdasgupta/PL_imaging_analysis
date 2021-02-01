from external_imports import *

led_curr = []
relative_intensity = []
with open("ledcal.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        led_curr.append(float(row[0]))
        relative_intensity.append(float(row[1]))
ledf = inter(led_curr,relative_intensity)

def parse(path):
    nominal_v_f = []
    led_current_f = []
    measured_current_f = []

    nominal_v_b = []
    led_current_b = []
    measured_current_b = []

    state = 'f'
    state_flipped = False
    with open(path+r"\iv.csv", 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if float(row[0]) <= nominal_v_f and state_flipped == False :
                state = 'b'
                state_flipped = True
            if state == 'f':
                nominal_v_f.append(float(row[0]))
                led_current_f.append(float(row[1]))
                measured_current_f.append(float(row[3]))
            else:
                nominal_v_b.append(float(row[0]))
                led_current_b.append(float(row[1]))
                measured_current_b.append(float(row[3]))
    return (nominal_v_f, led_current_f, measured_current_f), (nominal_v_b, led_current_b, measured_current_b)

def led_curr_at_1_sun(isc, LED_curr, meas_curr):
    f = inter(meas_curr, LED_curr)
    return f(isc)


def int_vs_curr(isc, area, LED_curr, meas_curr):
    led_curr_1sun = led_curr_at_1_sun(isc, LED_curr, meas_curr)
    intensities = ledf(LED_curr) / led_curr_1sun
    return intensities, np.array(meas_curr)/area

def reconstructed_jv(isc, voc, area,LED_curr, meas_curr):
    led_curr_1sun = led_curr_at_1_sun(isc, LED_curr, meas_curr)




