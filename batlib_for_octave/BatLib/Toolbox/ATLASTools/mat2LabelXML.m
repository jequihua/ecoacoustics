function [ output_args ] = mat2LabelXML( labels,name,externalChange,className,fileName )
%MAT2LABELXML Summary of this function goes here
%   Detailed explanation goes here

LabelTrack = com.mathworks.xml.XMLUtils.createDocument('LabelTrack');

root = LabelTrack.getDocumentElement;
root.setAttribute('classname',className);
root.setAttribute('externalchange',num2str(externalChange));
root.setAttribute('name',name);

for i=1:size(labels,2)
    label = LabelTrack.createElement('label');
    
    startTime=LabelTrack.createElement('starttime');
    startTime.appendChild(LabelTrack.createTextNode(num2str(labels(i).startTime)));
    label.appendChild(startTime);
    endTime=LabelTrack.createElement('endtime');
    endTime.appendChild(LabelTrack.createTextNode(num2str(labels(i).endTime)));
    label.appendChild(endTime);
    timestamp=LabelTrack.createElement('timestamp');
    timestamp.appendChild(LabelTrack.createTextNode(num2str(labels(i).timeStamp)));
    label.appendChild(timestamp);
    value=LabelTrack.createElement('value');
    value.appendChild(LabelTrack.createTextNode(num2str(labels(i).value)));
    label.appendChild(value);
    comment=LabelTrack.createElement('comment');
    comment.appendChild(LabelTrack.createTextNode(labels(i).comment));
    label.appendChild(comment);
    type=LabelTrack.createElement('type');
    type.appendChild(LabelTrack.createTextNode(labels(i).type));
    label.appendChild(type);
    text=LabelTrack.createElement('text');
    text.appendChild(LabelTrack.createTextNode(labels(i).text));
    label.appendChild(text);
    classentity=LabelTrack.createElement('classentity');
    classentity.appendChild(LabelTrack.createTextNode(labels(i).classentity));
    label.appendChild(classentity);
    
    root.appendChild(label);
    
end

xmlwrite(fileName,LabelTrack);

end

