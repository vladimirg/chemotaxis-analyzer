function angles = GetAnglesBetweenSteps(velocity_vectors)

[num_of_vecs, ~] = size(velocity_vectors);

angles = zeros(num_of_vecs - 1, 1);
for ix = 1:length(angles)
    vec1 = velocity_vectors(ix,:);
    vec2 = velocity_vectors(ix+1,:);
    angles(ix) = GetAngle(vec1, vec2);
end

end