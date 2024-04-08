import org.jocl.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import static org.jocl.CL.*;

// javac -classpath ".;lib\jocl-2.0.5.jar" KernalLauncher.java
// java -classpath ".;lib\jocl-2.0.5.jar" KernalLauncher

public class KernalLauncher {
    private static String readFile(String fileName) throws IOException {
        byte[] encoded = Files.readAllBytes(Paths.get(fileName));
        return new String(encoded, StandardCharsets.UTF_8);
    }

    public static void main(String[] args) throws IOException {
        CL.setExceptionsEnabled(true);
        final int width = 800;
        final int height = 800;
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        final int[] pixelData = new int[width * height];

        System.out.println("Rendering... (" + width + "px, " + height + "px)\n");
        String art = 
                "               O  o\n" + 
                "          _\\_   o\n" +
                "       \\\\/  o\\ .    this might take a while...\n" + 
                "       //\\___=\n" + 
                "          ''\n";
        System.out.println(art);

        // Obtain the number of platforms
        int numPlatformsArray[] = new int[1];
        CL.clGetPlatformIDs(0, null, numPlatformsArray);
        cl_platform_id platforms[] = new cl_platform_id[numPlatformsArray[0]];
        CL.clGetPlatformIDs(platforms.length, platforms, null);
        cl_platform_id platform = platforms[0];

        // Initialize the context properties
        cl_context_properties contextProperties = new cl_context_properties();
        contextProperties.addProperty(CL_CONTEXT_PLATFORM, platform);

        // Create a context for the GPU device
        cl_context context = CL.clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_GPU, null, null, null);

        // Get the list of GPU devices associated with the context
        long size[] = new long[1];
        CL.clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, null, size);
        cl_device_id devices[] = new cl_device_id[(int) (size[0] / Sizeof.cl_device_id)];
        CL.clGetContextInfo(context, CL_CONTEXT_DEVICES, size[0], Pointer.to(devices), null);

        // Create a command-queue
        cl_command_queue commandQueue = CL.clCreateCommandQueueWithProperties(context, devices[0], null, null);

        // Allocate the memory objects for the input- and output data
        cl_mem pixelBuffer = CL.clCreateBuffer(context, CL_MEM_WRITE_ONLY, Sizeof.cl_int * width * height, null, null);

        // Read the program source code and compile it
        String programSource = readFile("rayTracer.cl");
        cl_program program = CL.clCreateProgramWithSource(context, 1, new String[]{ programSource }, null, null);
        CL.clBuildProgram(program, 0, null, null, null, null);

        // Check for build errors
        long[] buildStatus = new long[1];
        CL.clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_STATUS, Sizeof.cl_long, Pointer.to(buildStatus), null);
        if (buildStatus[0] != CL_SUCCESS) {
            // If there's a build error, obtain and print the build log
            long[] logSize = new long[1];
            CL.clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0, null, logSize);
            byte[] logData = new byte[(int)logSize[0]];
            CL.clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, logSize[0], Pointer.to(logData), null);
            System.err.println("Build error:\n" + new String(logData, 0, logData.length - 1));
            System.exit(1);
        }

        // Create the kernel
        cl_kernel kernel = CL.clCreateKernel(program, "rayTracer", null);

        // Set the arguments for the kernel
        CL.clSetKernelArg(kernel, 0, Sizeof.cl_mem, Pointer.to(pixelBuffer));
        CL.clSetKernelArg(kernel, 1, Sizeof.cl_int, Pointer.to(new int[]{width}));
        CL.clSetKernelArg(kernel, 2, Sizeof.cl_int, Pointer.to(new int[]{height}));

        // Execute the kernel
        long global_work_size[] = new long[]{width, height};
        CL.clEnqueueNDRangeKernel(commandQueue, kernel, 2, null, global_work_size, null, 0, null, null);

        // Read the output data
        CL.clEnqueueReadBuffer(commandQueue, pixelBuffer, CL_TRUE, 0, width * height * Sizeof.cl_int, Pointer.to(pixelData), 0, null, null);

        image.setRGB(0, 0, width, height, pixelData, 0, width);

        // Save the image
        File outputfile = new File("image/image.png");
        ImageIO.write(image, "png", outputfile);
        System.out.println("Image successfully saved.");

        // Release kernel, program, and memory objects
        CL.clReleaseMemObject(pixelBuffer);
        CL.clReleaseKernel(kernel);
        CL.clReleaseProgram(program);
        CL.clReleaseCommandQueue(commandQueue);
        CL.clReleaseContext(context);
    }
}