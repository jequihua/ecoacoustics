function data = printout( data )

if data.printout.SonogramAudiofilesKey
    printMAINaudio
   
end
if data.printout.SonogramCallsKey
     printMAINcalls
end

if data.printout.WorksheetsKey
    compiletemplate(data)  
end


end

