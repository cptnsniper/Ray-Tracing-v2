inline int float3_color(float3 color) {
    // Convert each component to an integer in the range [0, 255]
    int ir = (int)(color.x * 255.0f);
    int ig = (int)(color.y * 255.0f);
    int ib = (int)(color.z * 255.0f);

    // Combine the components into a single int color in ARGB format
    // Assuming the alpha channel is fully opaque (255)
    int colorInt = (255 << 24) | (ir << 16) | (ig << 8) | ib;

    return colorInt;
}

__kernel void rayTracer(__global int* pixelData, const int width, const int height) {
    int i = get_global_id(0); // Column index
    int j = get_global_id(1); // Row index

    float u = (float)i / height * 2.0f - 1.0f;
    float v = (float)j / width * 2.0f - 1.0f;
    u *= width / height;
    v = 1.0f - v;

    if (i < width && j < height) { // Ensure to operate within image bounds
        float3 color = (float3)(u, v, 0.0f); // Define white color

        pixelData[j * width + i] = float3_color(color); // Set the pixel to white using the conversion function
    }
}