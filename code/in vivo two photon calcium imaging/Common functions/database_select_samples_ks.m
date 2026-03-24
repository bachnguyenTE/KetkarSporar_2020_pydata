% script to read in database and make various categories...

T=readtable('Summary.xlsx','ReadVariableNames',1);

stimbouts=T.stimbouts;
locnum=T.locnum;
inverted=T.invertedCellFlag;
stimcode=T.stimulusCode;
quality=T.qualitativeGoodness_0_10_;
driver=T.Driver;
moving=T.Moving;
layer=T.Layer_2_4_;
wavelength=T.StimulusWavelength;
responsivefly=T.responsivefly;

% reads now string for layers

% full field flashes
i_fff2s = strcmpi('LocalCircle_2sec_220deg_0degAz_0degEl_Sequential_LumDec_LumInc',(stimcode)); 
i_fff5s = strcmpi('LocalCircle_5sec_220deg_0degAz_0degEl_Sequential_LumDec_LumInc',(stimcode)); 
i_fff60s = strcmpi('LocalCircle_60sec_220deg_0degAz_0degEl_Sequential_LumDec_LumInc - Copy',(stimcode));

% ONgreyOFF stimuli
i_fff_01steps_5sONgrayOFF = strcmpi('LocalCircle_0.1steps_5sONgrayOFFgray_120deg_0degAz_0degEl_2Lum',(stimcode));

% 5 random steps
i_fff_05steps_10sONgrayOFF = strcmpi('LocalCircle_0.5steps_10s_randSteps_120deg_0degAz_0degEl_2Lum',(stimcode));

%Marvin's adaptation stimulus
i_FullField_OFF_2s_OFF_2s = strcmpi('FullField_OFF_2s_OFF_2s_1s_BG_15s_con_m25_lum_0d75to0_NonRand',(stimcode));
% i_FullField_ON_3s_OFF_3s = strcmpi('FullField_OFF_2s_OFF_2s_1s_BG_15s_con_m25_lum_0d75to0_NonRand',(stimcode));

% Slow ramp down 100 steps
i_slow_ramp_down_500ms_100 = strcmpi('slow_ramp_down_500ms_100',(stimcode));
i_slow_ramp_down_6s_100 = strcmpi('slow_ramp_down_6s_100',(stimcode));
i_LocalCircle_5min_control_for_the_ramp = strcmpi('LocalCircle_5min_control_for_the_ramp',(stimcode));

% 30 s adapting step, then OFF step to the same luminance
i_Diff_Adapting_OFFstep_6s = strcmpi('Diff_Adapting_OFFstep_6s',(stimcode));


i_15bouts = stimbouts>=15;
i_inverted = (inverted == 1);
i_moving = (moving == 1);
i_active = (quality ~= 0);
i_responsivefly = (responsivefly == 1);  %ms addition

%changed to read strings
% i_L4M2 = strcmpi('L4M2',layer);
% i_L4M4 = strcmpi('L4M4',layer);
i_lob1 = strcmpi('lob1',layer);
i_447 = (wavelength == 447);


i_10to15min = locnum>10&locnum<15;
i_15to20min = locnum>15&locnum<20;
i_20to25min = locnum>20&locnum<25;
i_25to30min = locnum>25&locnum<30;
i_30to35min = locnum>30&locnum<35;

i_20to30min = locnum>20&locnum<30;
i_recovery = locnum>60;


% L2 / L3 >>UAS-GCaMP6f
i_UASGCaMP6F_L221DhhGal4 = strcmpi('UASGCaMP6F_L221DhhGal4',driver);
i_UASGCaMP6F_L3MH56Gal4 = strcmpi('UASGCaMP6F_L3MH56Gal4',driver);

%UAS-GCaMP6f;L2[21Dhh-Gal4] or MH56-Gal4 x w+
i_UASGCaMP6f_L3MH56Gal4_cross_to_w = strcmpi('UASGCaMP6f_L3MH56Gal4_cross_to_w',driver);
i_UASGCaMP6F_L221DhhGal4_cross_to_w = strcmpi('UASGCaMP6F_L221DhhGal4_cross_to_w',driver);

i_UASGCaMP6f_L3MH56Gal4_cross_to_w_mars = strcmpi('UASGCaMP6f_L3MH56Gal4_cross_to_w_mars',driver);
i_UASGCaMP6F_L221DhhGal4_cross_to_w_mars = strcmpi('UASGCaMP6F_L221DhhGal4_cross_to_w_mars',driver);

% Voltage sensors ASAP2
i_L3MH56Gal4_Asap2f = strcmpi('L3MH56Gal4_Asap2f',driver);
i_L221DhhGal4_Asap2f = strcmpi('L221DhhGal4_Asap2f',driver);

%current 150 - for the ND filters experiment
i_UASGCaMP6f_MH56Gal4_L3_150_ND06 = strcmpi('UASGCaMP6f_L221DhhGal4_MH56Gal4_L3_150_ND06',driver);
i_UASGCaMP6f_MH56Gal4_L3_150_ND13 = strcmpi('UASGCaMP6f_L221DhhGal4_MH56Gal4_L3_150_ND13',driver);
