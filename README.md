# rc2atomneb

`rc2atomneb.py` converts recombination-line source data into
AtomNeb-style FITS files for the `atomic-data-rc/` directory of [AtomNeb-py](https://github.com/atomneb/AtomNeb-py).

It is a Python analogue of the IDL generator `gen_rc_atomneb.pro` used by [AtomNeb-idl](https://github.com/atomneb/AtomNeb-idl).

## Requirements

```bash
pip install numpy astropy
```

## Input

Expected source tree:

```text
rc_data/
  rc_collection/
  rc_SH95/
  rc_PPB91/
  rc_he_ii_PFSD12/
  rc_n_iii_FSL13/
  rc_o_iii_SSB17/
```

The attached short `rc_o_iii_SSB17/OIIlines_ABC.txt` is enough for parser tests,
but it is not the complete SSB17 dataset. To reproduce the full AtomNeb SSB17
FITS files, replace it with the full CDS file `OIIlines_ABC` from VizieR VI/150.

The recombination data are as follows:

* *rc_collection*, effective recombination coefficients for C II ([Davey et al. 2000](http://adsabs.harvard.edu/abs/2000A%26AS..142...85D)), N II ([Escalante and Victor 1990](http://adsabs.harvard.edu/abs/1990ApJS...73..513E)), O II ([Storey 1994](http://adsabs.harvard.edu/abs/1994A%26A...282..999S); [Liu et al. 1995](http://adsabs.harvard.edu/abs/1995MNRAS.272..369L)), and Ne II ions ([Kisielius et al. 1998](http://adsabs.harvard.edu/abs/1998A%26AS..133..257K)), including Branching ratios (Br) for O II and N II ions.

* *rc_SH95*, hydrogenic ions for Z=1 to 8, namely H I, He II, Li III, Be IV, B V, C VI, N VII, and O VIII ions from [Storey and Hummer (1995)](http://adsabs.harvard.edu/abs/1995MNRAS.272...41S).

* *rc_PPB91*, effective recombination coefficients for H, He, C, N, O, Ne ions from [Pequignot, Petitjean and Boisson (1991)](http://adsabs.harvard.edu/abs/1991A%26A...251..680P).

* *rc_he_ii_PFSD12*, effective He I recombination coefficients from [Porter et al (2012)](
http://adsabs.harvard.edu/abs/2012MNRAS.425L..28P) and [(2013)](http://adsabs.harvard.edu/abs/2013MNRAS.433L..89P).

* *rc_n_iii_FSL13*, effective N II recombination coefficients from [Fang, Storey and Liu (2011)](
http://adsabs.harvard.edu/abs/2011A%26A...530A..18F) and [(2013)](http://adsabs.harvard.edu/abs/2013A%26A...550C...2F).

* *rc_o_iii_SSB17*, effective O II recombination coefficients of 8889 recombination lines for Cases A, B, and C from [Storey, Sochi and Bastin (2017)](
http://adsabs.harvard.edu/abs/2017MNRAS.470..379S).


## Output

By default:

```bash
python rc2atomneb.py --rc-root ./rc_data --out-dir ./atomic-data-rc
```

writes:

```text
atomic-data-rc/
  rc_collection.fits
  rc_PPB91.fits
  rc_SH95.fits
  rc_he_ii_PFSD12.fits
  rc_n_iii_FSL13.fits
  rc_o_iii_SSB17.fits
  rc_o_iii_SSB17_orl_case_b.fits
```

## FITS layout

Each FITS file follows the recombination-data layout used by
`gen_rc_atomneb.pro` in [AtomNeb-idl](https://github.com/atomneb/AtomNeb-idl):

```text
PRIMARY
List
References
data extensions...
```

`List` contains at least:

```text
Aeff_Data
Extension
```

and some collections include additional line metadata. `References` contains:

```text
AtomicData
Reference
```

Data extensions are named from the `Aeff_Data` field, e.g.

```text
c_iii_aeff
h_ii_aeff_a
he_ii_aeff_PFSD12
n_iii_aeff_1
o_iii_aeff_b_1234
```

## Build selected collections

```bash
python rc2atomneb_v3_optional_ssb17.py \
  --rc-root ./rc_data \
  --out-dir ./atomic-data-rc \
  --collections collection,ppb91,sh95
```

Accepted collection names:

```text
all
collection
ppb91
sh95
pfsd12
fsl13
ssb17
```

## SSB17 options

Make only the compact optical Case B file:

```bash
python rc2atomneb_v3_optional_ssb17.py \
  --rc-root ./rc_data \
  --out-dir ./atomic-data-rc \
  --collections ssb17 \
  --ssb17-no-full
```

Make only the full SSB17 file:

```bash
python rc2atomneb_v3_optional_ssb17.py \
  --rc-root ./rc_data \
  --out-dir ./atomic-data-rc \
  --collections ssb17 \
  --ssb17-no-case-b-optical
```

## Notes

- The script writes FITS image extensions for grid-like effective recombination
  coefficient arrays, following the IDL `mwrfits` behavior.
- The script writes FITS binary table extensions for line lists, references,
  analytic-fit coefficients, wavelength tables, and branching-ratio tables.
- The SSB17 parser supports the full CDS `OIIlines_ABC` file structure and can
  also parse shortened test files, but shortened files cannot reproduce the full
  AtomNeb SSB17 products.


## Optional SSB17 O II products

The SSB17 O II recombination products are **not created by default** because the
main source table is very large and may not be included in small test archives.

By default:

```bash
python rc2atomneb_v3_optional_ssb17.py \
  --rc-root ./rc_data \
  --out-dir ./atomic-data-rc
```

or:

```bash
python rc2atomneb_v3_optional_ssb17.py \
  --rc-root ./rc_data \
  --out-dir ./atomic-data-rc \
  --collections all
```

creates only:

```text
rc_collection.fits
rc_PPB91.fits
rc_SH95.fits
rc_he_ii_PFSD12.fits
rc_n_iii_FSL13.fits
```

It does **not** create:

```text
rc_o_iii_SSB17.fits
rc_o_iii_SSB17_orl_case_b.fits
```

To request the SSB17 products explicitly:

```bash
python rc2atomneb_v3_optional_ssb17.py \
  --rc-root ./rc_data \
  --out-dir ./atomic-data-rc \
  --collections ssb17
```

or together with the default products:

```bash
python rc2atomneb_v3_optional_ssb17.py \
  --rc-root ./rc_data \
  --out-dir ./atomic-data-rc \
  --collections all,ssb17
```

When SSB17 is requested, the script checks for one of these files:

```text
rc_data/rc_o_iii_SSB17/OIIlines_ABC.txt
rc_data/rc_o_iii_SSB17/OIIlines_ABC
```

If neither file exists, download the full table from CDS:

```text
https://cdsarc.cds.unistra.fr/ftp/VI/150/DataFiles/OIIlines_ABC
```

or from the catalog page:

```text
https://cdsarc.cds.unistra.fr/viz-bin/cat/VI/150
```

Then place it in:

```text
rc_data/rc_o_iii_SSB17/OIIlines_ABC
```

or:

```text
rc_data/rc_o_iii_SSB17/OIIlines_ABC.txt
```

The shortened `OIIlines_ABC.txt` used for parser testing is not sufficient to
reproduce the complete AtomNeb SSB17 FITS products.
