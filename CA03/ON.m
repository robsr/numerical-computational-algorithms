function Q=ON(A)

sz=size(A);
n=sz(1);
q=cell(1,n);

q{1}=A(:,1)/norm(A(:,1));
Q=q{1};
for i=2:n
    q{i}= A(:,i);
    for j=1:i-1
        q{i} = q{i}-( A(:,i)'*(q{j}) )*q{j};
    end
    q{i}=q{i}/norm(q{i});
    Q=[Q q{i}];
end
end
