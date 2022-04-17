function w = getWeights(ref)
M = containers.Map();
amino_acids = geneticcode();
codons = revgeneticcode();
f1 = fieldnames(codons);
f2 = fieldnames(amino_acids);

for i = 1:length(f2) %finding the synonymous codons
    amino_acid = amino_acids.(f2{i});
    try
        M(f2{i}) = codons.(amino_acid);
    catch
    
    end
end


counts = codoncount(ref);
f3 = fieldnames(counts);

for i = 1:length(f3)
    if counts.(f3{i}) == 0
        counts.(f3{i}) = 0.5;
    end

end

f4 = keys(M);
frequent = containers.Map();

for i = 1:length(f4) %calculating the frequencies
    k = f4{i};
    sumF = 0;
    sCodons = M(k);
    for j = 1:length(sCodons)
        sumF = sumF + counts.(sCodons{j});
    end
     w1 = length(sCodons) .^ -1 * sumF;
     n1 = counts.(k);
    frequent(k) = n1/w1;
end
disp(frequent)

w = containers.Map();

for i = 1:length(f4) %calculating the weights
    k = f4{i};
    maxCodon = 0.5;
    sCodons = M(k);
    for j = 1:length(sCodons)
        if frequent(sCodons{j}) > maxCodon
            maxCodon = frequent(sCodons{j});
        end
    end
    w(k) = frequent(k) / maxCodon;
end