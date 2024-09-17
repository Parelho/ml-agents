import time
import space_state_legs as cp
import json

freq_hz = 1000
periodo = 1 / freq_hz
pitch = 0.8
ad = [ 0.5 , 1.0 , -0.5 , 0.3]
pe = [ 0.0 , 0.0]

def executa_perna():
    global periodo, pitch, ad, pe
    torques = cp.executa_perna(pitch, ad[0], ad[2], ad[1], ad[3], pe[0], pe[1], periodo)
    return {
        "t0":torques[0],
        "t1":torques[1],
        "t2":torques[2],
        "t3":torques[3]
    }


if __name__ == "__main__":
    result = executa_perna()
    with open("result.json", "w") as f:
        json.dump(result, f)