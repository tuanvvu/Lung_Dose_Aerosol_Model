## Calculating the lung dose of aerosols in different micro-environments

### Deposition fraction of particle in the human respiratory system
This work provides the regional lung deposition calculation program of aeerosols based on the recommended deposition formulae and their parametric values derived from the International Commission on Radiological Protection ( ICRP, Publication 66). The ICRP model is a semi-empirical model which determines the deposited fraction of particles in five regions of the airway system (the nose and mouth, throat and larynx, upper airways, lower airways, and alveoli) using both numerical fitting of experimental data and theoretical calculations. The deposition of particles is controlled by different transport processes, which strongly depend on particle size. In the empirical modelling of the deposition data, the deposited fraction of particles is controlled by two kinds of deposition processes known as aerodynamic and thermodynamic transport: thermodynamic transport predominantly accounts for the deposition of small particles (<0.1 µm), while aerodynamic transport mainly influences large particles (>1.0 µm). For particles in a transition regime size range (0.1-1.0 µm), deposition efficiencies are determined by both aerodynamic and thermodynamic transport processes.

Input properties of particles into ICRP model:
1. Particle's properties:
- Particle size (in nm, model valid from a wide size range of 2 nm- 20 um)
- Particle hygroscopic's growth factor
- Particle effective density & shape (dimension) factor
2. Subject's lung structure parameters/ Human activites: FR, respiration frequency; VE, minute ventilation; VT, tidal volume of exposed subject; Fn, fraction of total ventilator airflow passing through the nose; FRC, functional residual capacity of the exposed subject; VdET, anatomical dead space of the ET; VdBB, anatomical dead space of the trachea and bronchi; Vdbb, anatomical dead space of the bronchioles; SFt, SFb, SFa, scaling parameters for males; V, volumetric flow rate of inspired air.

### Lung Dose Calculation
The regional lung dose is defined as the proportion of inhaled particles deposited in the respiratory tract during an exposure time period (Δt) and is calculated from the following:
Dosei=∑VE×C×DFi×Δt
where DFi is the total deposited fraction of particles in the lung region (i). C is the total concentration of particles by number, surface area or mass. VE is defined as the ventilation rate which is the volume of gas inhaled or exhaled from the lungs during a unit time period (m3/h), and Δt is the exposure time duration.

### Refererences:
1. Vu, T.V., Zauli-Sajani, S., Poluzzi, V. et al. Air Qual Atmos Health (2018) 11: 615. Factors controlling the lung dose of road traffic-generated sub-micrometre aerosols from outdoor to indoor environments. https://doi.org/10.1007/s11869-018-0568-2
2. Vu, T.V., Ondracek, J., Zdímal, V. et al. Air Qual Atmos Health (2017) 10: 1. Physical properties and lung deposition of particles emitted from five major indoor sources. https://doi.org/10.1007/s11869-016-0424-1
3. Vu, T.V., Delgado-Saborit, J.M. & Harrison, R.M. Air Qual Atmos Health (2015) 8: 429.A review of hygroscopic growth factors of submicron aerosols from different sources and its implication for calculation of lung deposition efficiency of ambient aerosols. https://doi.org/10.1007/s11869-015-0365-0

