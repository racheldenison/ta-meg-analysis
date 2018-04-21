% donders_playground_sourcerecon.m

%% setup
mrifile = '/Volumes/RACHO/Hypatia_20160416/Anatomy/mr19900405/mr/mri/T1.mgz';
headshapefile = '/Volumes/RACHO/Hypatia_20160416/TADetectDiscrim/MEG/R0983_20150813/R0983_8.13.15_HS.txt';
% headshapedir = 

%% load mri
mri = ft_read_mri(mrifile, 'dataformat', 'freesurfer_mgz');

%% realign
cfg = [];
cfg.method = 'headshape';
cfg.headshape.headshape = headshapefile;
mrir = ft_volumerealign(mri);