function mpc = case95
% CASE9 + CASE5    Power flow data for 9 bus, 3 generator case. + case5
% MATPOWER case9 + case5 data
% http://www.pserc.cornell.edu//matpower/

% R. D. Zimmerman, C. E. Murillo-Sánchez, and R. J. Thomas,
% "MATPOWER: Steady-State Operations, Planning and Analysis Tools
% for Power Systems Research and Education," Power Systems, IEEE
% Transactions on, vol. 26, no. 1, pp. 12-19, Feb. 2011

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	345	1	1.1	0.9;
	2	2	0	0	0	0	1	1	0	345	1	1.1	0.9;
	3	2	0	0	0	0	1	1	0	345	1	1.1	0.9;
	4	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
	5	1	90	30	0	0	1	1	0	345	1	1.1	0.9;
	6	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
	7	1	100	35	0	0	1	1	0	345	1	1.1	0.9;
	8	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
	9	1	125	50	0	0	1	1	0	345	1	1.1	0.9;
    10	2	0	0	0	0	1	1	0	230	1	1.1	0.9;
	11	1	300	98.61	0	0	1	1	0	230	1	1.1	0.9;
	12	2	300	98.61	0	0	1	1	0	230	1	1.1	0.9;
	13	1	400	131.47	0	0	1	1	0	230	1	1.1	0.9;
	14	2	0	0	0	0	1	1	0	230	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	300	-300	1	100	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
	2	163	0	300	-300	1	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
	3	85	0	300	-300	1	100	1	270	10	0	0	0	0	0	0	0	0	0	0	0;
    10	40	0	30	-30	1	100	1	40	0	0	0	0	0	0	0	0	0	0	0	0;
	10	170	0	127.5	-127.5	1	100	1	170	0	0	0	0	0	0	0	0	0	0	0	0;
	12	323.49	0	390	-390	1	100	1	520	0	0	0	0	0	0	0	0	0	0	0	0;
	13	0	0	150	-150	1	100	1	200	0	0	0	0	0	0	0	0	0	0	0	0;
	14	466.51	0	450	-450	1	100	1	600	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	4	0	0.0576	0	250	250	250	0	0	1	-360	360;
	4	5	0.017	0.092	0.158	250	250	250	0	0	1	-360	360;
	5	6	0.039	0.17	0.358	150	150	150	0	0	1	-360	360;
	3	6	0	0.0586	0	300	300	300	0	0	1	-360	360;
	6	7	0.0119	0.1008	0.209	150	150	150	0	0	1	-360	360;
	7	8	0.0085	0.072	0.149	250	250	250	0	0	1	-360	360;
	8	2	0	0.0625	0	250	250	250	0	0	1	-360	360;
	8	9	0.032	0.161	0.306	250	250	250	0	0	1	-360	360;
	9	4	0.01	0.085	0.176	250	250	250	0	0	1	-360	360;
  10	11	0.00281	0.0281	0.00712	400	400	400	0	0	1	-360	360;
  10	13	0.00304	0.0304	0.00658	0	0	0	0	0	1	-360	360;
  10	14	0.00064	0.0064	0.03126	0	0	0	0	0	1	-360	360;
  11	12	0.00108	0.0108	0.01852	0	0	0	0	0	1	-360	360;
  12	13	0.00297	0.0297	0.00674	0	0	0	0	0	1	-360	360;
  13	14	0.00297	0.0297	0.00674	240	240	240	0	0	1	-360	360;
	7   14   0.01	0.085	0.176	250	250	250	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% area data
%	area	refbus
mpc.areas = [
	1	5;
];

%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	1500	0	3	0.11	5	150;
	2	2000	0	3	0.085	1.2	600;
	2	3000	0	3	0.1225	1	335;
];
