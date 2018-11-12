# Nomenclature
* Event: Event that releases chemical element (e.g., supernova).
* Progenitor: Progenitor object that will eventually create an event (e.g., massive star).
* Delta: Time interval between two consecutive events.
* Gamma: Time interval between the formation of two consecutive progenitors.

### List of event times

- - - - - 
**tEvents_PDF\_???\_???\_gamma\_?????\_tend_\????.in**

- PDF: Type of probability distribution function with time boundaries [tmin, tmax].
	- box: Square box where each time can be selected with equal probability.
	- pow: Power law in the form of 1/t.
- gamma: Constant time interval between the formation of progenitors.
- tend: Formation time of the last progenitor.
	- **Note**: This is not the final time value in the file. The last time value is associated with the last enriching event, which is likely to occur after tend.
	
Example: tEvents_pow_1e7_1e10_gamma_1_00e6_tend_15e9.in

- Power-law PDF sampled from 1e7 to 1e10 yr, gamma = 1e6 yr, last progenitor formed at 15e9 yr.
	
- - - - - 