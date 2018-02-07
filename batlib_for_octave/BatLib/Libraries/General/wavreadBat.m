function values = wavreadBat(fileName)
% function values = wavreadBat(fileName)
% This function reads common .wav file formats used in bat experiments.
% Input:  fileName--character array of name of one .wav file, including '.wav'
% Output: values.x--Nx1 vector of audio data (or Nx2 for stereo, Nx4 for 4-ch data).
%                   8-bit: 0..255, returned as uint8 for memory limitations (subtract 128 when computing).
%                   16-bit: -32768..32767, returned as int16.
%         values.fs--sampling rate, Hz (may or may not be time-expanded fs)
%         values.source--Estimated source of bat recording (Pettersson, Avisoft, Titley, Sonobat, Batscan, or unknown)
%         values.bitsPerSample--number of bits per audio sample
%         values.noiseFloordB--estimate of noise floor for int8 or int16 integer data, in dB
%         values.peakSNRdB--estimate of signal-to-noise ratio, in dB: peak energy - noiseFloordB
%         values.audioFormat--1==PCM, 2==Microsoft ADPCM.  Only PCM is currently supported (cannot decode Batscan files).
%                             See http://www.iana.org/assignments/wave-avi-codec-registry for complete list.
%         values.numChannels--number of audio streams in data (1,2, or 4 only)
%         values.riff, .wave, .fmt, .data--sanity check of text present in all .wav headers.  If
%                                          text doesn't match field name, something is wrong...
%         values.lastRIFF--byte number of 'R' in 'RIFF' in header.  Some sources (e.g., Pettersson) have 3 or 4 'RIFF's present.
%         values.chunkSize--total file length (starting from 'R' in last 'RIFF') minus 8 bytes
%         values.subchunk1Size--Number of bytes between 'fmt ' and 'data' text tags (usually 16)
%         values.subchunk2Size--Number of bytes in audio data
%         values.dataPtr: Byte number of 'd' in 'data' tag (offset from lastRIFF value)
% For details on .wav header, see the following: http://ccrma.stanford.edu/CCRMA/Courses/422/projects/WaveFormat/
% Note: values.x = [] for Batscan files because the MS ADPCM codec is not currently supported.

% Mark Skowronski, July 6, 2005
% Comments or bug fixes?  Send to markskow@cnel.ufl.edu

% Check file for validity:
fid = fopen(fileName);
if fid==-1,
   error('ERROR: invalid file name.  Quitting.');
end;

% Read header:
fseek(fid,0,-1); % Move cursor to beginning of file
head = fread(fid,5000,'uint8'); % Read as unsigned bytes

% Find last RIFF:
h = char(head(:)'); % Convert to ASCII text
r = findstr('RIFF',h); % At least one, always
if isempty(r),
   error('ERROR: not a valid RIFF file.  Quitting.');
end;
r = r(end); % Point to last one (normal .wav files and Avisoft only have one, but Pettersson have 3 or 4)

% Point at last RIFF, header follows:
head = head(r:end);

% Pull out various parameters:
values = struct([]);
values(1).lastRIFF = r;
values.riff = char(head(1:4))'; % 'RIFF' in ASCII, always
values.chunkSize = head([1:4]+4); % Chunk size, in little endian
values.chunkSize = 16.^[0:2:6]*values.chunkSize(:); % Convert 4-byte little endian to number
values.wave = char(head([1:4]+8))'; % 'WAVE' in ASCII
values.fmt = char(head([1:4]+12))'; % 'fmt ' in ASCII
values.subchunk1Size = head([1:4]+16); % 16 bytes for PCM, 18 bytes for Pettersson (2 extra after BitsPerSample), sometimes larger when extra parameters are defined between bitsPerSample and subChunk2ID
values.subchunk1Size = 16.^[0:2:6]*values.subchunk1Size(:); % Convert 4-byte little endian to number
values.audioFormat = head([1:2]+20);
values.audioFormat = 16.^[0,2]*values.audioFormat(:); % Convert 2-byte little endian to number
values.numChannels = head([1:2]+22);
values.numChannels = 16.^[0,2]*values.numChannels(:);
values.fs = head([1:4]+24);
values.fs = 16.^[0:2:6]*values.fs(:);
values.bitsPerSample = head([1:2]+34);
values.bitsPerSample = 16.^[0,2]*values.bitsPerSample(:);

% Find 'data' field (varies since there could be extra parameters present after bitsPerSample):
dataPtr = findstr('data',char(head(:)')); % one, always (unless also used in a comment field)
dataPtr = dataPtr(end); % Use last one (could be 'data' in some comments before 'data' field
values.dataPtr = dataPtr;
values.data = char(head([1:4]+values.dataPtr-1))';
values.subchunk2Size = head([1:4]+values.dataPtr-1+4);
values.subchunk2Size = 16.^[0:2:6]*values.subchunk2Size(:); % Convert 4-byte little endian to number

if values.audioFormat==1, % PCM (no codec)
   % Get data (8- and 16-bit data treated differently, according to RIFF specification):
   fseek(fid,values.dataPtr-1+values.lastRIFF-1,-1); % Move cursor to beginning of data
   if values.bitsPerSample==8, % Read as 8-bit unsigned integers, return as uint8 variable, not double precision
      values.x = fread(fid,values.subchunk2Size,'*uint8');
   elseif values.bitsPerSample==16, % Read as 16-bit signed integers, return as int16 variable, not double precision
      values.x = fread(fid,values.subchunk2Size/2,'*int16'); % 16-bit signed integers, 2's complement
   end;

   % Write interlaced data as multiple channels, if needed:
   if values.numChannels==2,
      values.x = [values.x(1:2:end),values.x(2:2:end)]; % Left, right channels
   elseif values.numChannels==4,
      values.x = [values.x(1:4:end),values.x(2:4:end),values.x(3:4:end),values.x(4:4:end)];
   end;

   % Estimate noise floor:
   fs = values.fs;
   if fs<100000,
      fs = fs*10; % Time expanded, so correct
   end;
   numFrames = floor(length(values.x)/fs/(2/1000)); % 2 ms windows to estimate power
   g = zeros(1,numFrames);
   for p=1:numFrames,
      xIndex = [1:floor(fs*2/1000)]+(p-1)*floor(fs*2/1000); % 2 ms frames, no overlap
      if values.bitsPerSample==8,
         g(p) = mean((double(values.x(xIndex))-128).^2); % mean squared value
      else
         g(p) = mean(double(values.x(xIndex)).^2); % mean squared value
      end;
   end;
   g = 10*log10(g+eps); % rms, in dB (eps added in case g(p)=0 for some p)
   [aa,bb] = hist(g,ceil(numFrames/10));
   [junk,bbIndex] = max(aa);
   values.noiseFloordB = bb(bbIndex); % dB estimate
   values.peakSNRdB = max(g)-bb(bbIndex);
else % No codec support at this time
   values.x = []; % No data stored
   values.noiseFloordB = []; % No noise floor estimate
   values.peakSNRdB = []; % No peak SNR estimate
end;

% Determine source:
fseek(fid,-100,1); % Move cursor 100 bytes from end of file
foot = fread(fid,inf,'uint8'); % Read as unsigned bytes
if ~isempty(findstr('Pettersson',h(:)')), % Pettersson
   values.source = 'Pettersson';
elseif ~isempty(findstr('MMMMMMMMM',h(:)')), % Sonobat
   values.source = 'Sonobat';
elseif ~isempty(findstr('fact',h(:)')) & values.audioFormat==1, % Titley
   values.source = 'Titley'; % 'fact' is normally found when audioFormat>1
elseif ~isempty(findstr('TIME',foot(:)')),
   values.source = 'Avisoft';
elseif values.audioFormat==2, % ADPCM codec, probably from Batbox/Batscan
   values.source = 'Batscan';
else
   values.source = 'Unknown source';
end;

% Close:
fclose(fid);




% Bye!