# pnp-aoa
PNP solver for studying the transport in aoa

## Usage

```
pnp_num( charge_strgengh, enzymatic_rate, bulk_concentration );
```
The simulation result will be saved in a mat file.

## Example

```
pnp_num( 1, 5, 10 );
```

for charges on, rate 5 per second, and bulk concentration 10 nM.

Note: though the input of the bluk concentration is in the form of moles/L, the concentration profiles in the output file are of moles/m^3.
