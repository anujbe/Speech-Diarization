function [out]=bic_func(str1,end1,str2,end2,av,lambda)

sz1=end1-str1+1;
sz2=end2-str2+1;
sz3=sz1+sz2;

z1=zeros(sz1,13);
z2=zeros(sz2,13);
z3=zeros(sz3,13);

z1_cov=zeros(13,13);
z2_cov=zeros(13,13);
z3_cov=zeros(13,13);

for i=str1:end1
    z1(i-str1+1,:)=av(i,:);
end

for i=str2:end2
    z2(i-str2+1,:)=av(i,:);
end

for i=1:sz3
    if(i>sz1)
        z3(i,:)=av((str2-1)+(i-sz1),:);
    else
        z3(i,:)=av((str1-1)+i,:);
    end
end

z1_cov(:,:)=cov(z1);
z2_cov(:,:)=cov(z2);
z3_cov(:,:)=cov(z3);


 C=log(det(z3_cov));
 B=log(det(z2_cov));
 A=log(det(z1_cov));
 R=(sz3.*C) - (sz2.*B) - (sz1.*A);
 P=.5.*(13 + 0.5.*13.*(13+1)).*log(sz3);
 out = - R + lambda.*P ;
