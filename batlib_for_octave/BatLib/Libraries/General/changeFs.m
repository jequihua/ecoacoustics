function data = changeFs(data,oldFs,newFs)

i = data.audioinfo.SamplingFrequency==oldFs;
data.audioinfo.SamplingFrequency(i) = newFs;
end