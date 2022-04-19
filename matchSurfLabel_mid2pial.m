function [label_ds,pialv_ds,midv_ds_updated,midx]=matchSurfLabel_mid2pial(label,pialv,midv,midv_ds)
%use this function over matchSurfLabel if you want to relabel a downsampled
%pial surface. Problem with using matchSurfLabel is that sup temp and inf
%frontal might get mislabeled as they are too close together.

if size(midv,1)~= size(pialv,1)
    error('pial and mid surface do not appear to be matched!')
end

np=size(midv_ds,1);
label_ds=zeros(np,1);
midv_ds_updated=midv_ds;
pialv_ds=zeros(size(midv_ds));
midx=zeros(size(midv_ds,1),1);

steps=1000;%how many points we can handle in parallel (to speed up the loop)

for i=1:steps:np
    if i+steps-1>np
        ei=np;
    else
        ei=i+steps-1;
    end
    
%     disp([num2str(i) ' of ' num2str(np)])
    cp=midv_ds(i:ei,:);
    
    d1=midv(:,1)-cp(:,1)';
    d2=midv(:,2)-cp(:,2)';
    d3=midv(:,3)-cp(:,3)';
    
    d=sqrt(d1.^2+d2.^2+d3.^2);
    
    [~,midx(i:ei)]=min(d);

    label_ds(i:ei)=label(midx(i:ei));
    midv_ds_updated(i:ei,:)=midv(midx(i:ei),:);
    pialv_ds(i:ei,:)=pialv(midx(i:ei),:);%uses the fact that pial and mid surface points are 1:1 matched
end