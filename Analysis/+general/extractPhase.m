% Extract instantaneous phase from a signal
%
% This function extracts instantaneous phase of a signal based on Hilbert transform.
%
%  USAGE
%   phase_deg = general.extractPhase(eeg)
%   eeg             Vector that contains a signal
%   phase_deg       Vector that contains extracted phase. Vector values are in degrees.
%
function [phase_deg] = extractPhase(eeg)
%     phase_rad = phase(hilbert(eeg));
    phase_rad = angle(hilbert(eeg)); % faster than phase
    phase_deg = rad2deg(mod(phase_rad, 2*pi));
end