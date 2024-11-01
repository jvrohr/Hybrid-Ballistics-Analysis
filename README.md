# Hybrid-Ballistics-Analysis

This is a simple simulator for the internal ballistics of a hybrid rocket engine/motor with a feeding system based on blowdown of auto-pressurized nitrous oxide (N2O) as oxidizer and paraffin (C73H124) as fuel. Nitrous oxide properties were taken from ESDU 91022. User has two options of runs:

Blowdown: evacuation of nitrous oxide from the tank to open air, so the downstream pressure is always equal to the atmospheric pressure. 
Burn: a normal burn accounting for the combustion chamber pressure in the blowdown model.

Feel free to modify it and use your own way. An example file is provided with how to run a simple motor blowdown and burn with plots (```example.py```).

## Possible improvements

- Add more fuel and oxidizer options.
- Add a nozzle model to get thrust, etc. By now it just gets the Cf from CEA and calculates the thrust.
- Better account for transients as for amateur engines they are a considerable amount of time.
- Model the feeding system lines pressure drop.
