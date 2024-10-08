# Dependencies installation

- [1. Install dependencies: `Docker`](#1-install-dependencies-docker)
   * [1.1. Install Docker](#11-install-docker)
   * [1.2. Pull the Docker image](#12-pull-the-docker-image)
   * [1.3. Create a Docker container](#13-create-a-docker-container)
   * [1.4. Use the Docker container](#14-use-the-docker-container)
   * [1.5. Editing files](#15-editing-files)
- [2. Test the installation](#2-test-the-installation)

---

# 1. Install dependencies: `Docker`

## 1.1. Install Docker
As a first step, install [Docker Desktop](https://www.docker.com/products/docker-desktop/). You can follow the instruction on the [official guide](https://docs.docker.com/get-docker/) and the [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/). Please read the guides thoroughly.

## 1.2. Pull the Docker image
From a terminal, run the command

```bash
docker pull pcafrica/adv_sci_comp
```

You can check your images with `docker image ls`.

## 1.3. Create a Docker container
To use your image you need to create a Docker container out of it. The following command creates a container with the image just downloaded, gives it an easy to remember name (`--name adv_sci_comp`) and share a folder with the host so that files can be transferred easily (`-v /path/to/host/folder:/home/jellyfish/shared-folder`):

```bash
docker run --name adv_sci_comp -v /path/to/host/folder:/home/dealii/shared-folder -it -d pcafrica/adv_sci_comp
```

> **Note**: to avoid problems, paths should not contain white spaces or special characters.

You have now created a container. To turn on the container type:

```bash
docker start adv_sci_comp
```

## 1.4. Use the Docker container
To enter into the container run:

```bash
docker exec -it adv_sci_comp /bin/bash
```

You can leave the container and return to your OS with `exit`. You can check your containers and their status with the command

```bash
docker ps -a
```

If the status of the container is `Up`, you can stop it with

```bash
docker stop adv_sci_comp
```

Once you have created your container remember to **do not** use again the commad `run` but just `start`. Otherwise you will create every time a new container. If you want to remove a container you created for mistake you can run:

```bash
docker rm <name-of-the-container>
```

> **Note**: if the Docker container returns "Permission denied." errors when trying to create or modify files, follow these instructions (only once):
>    ```bash
>    # On the host.
>    id -u $(whoami) # Get your user ID.
>    id -g $(whoami) # Get your group ID.
>
>    docker exec --user root -it adv_sci_comp /bin/bash
>    /sbin/usermod -u 1000 dealii # Replace 1000 with your user ID.
>    /sbin/groupmod -g 1000 dealii # Replace 1000 with your group ID.
>    exit
>    ```
>    After that, you should be able to execute the container normally.

## 1.5. Editing files

1. You only need a text editor or IDE installed on the host OS (or a non-graphical editor or editor inside Docker) of your choice. [Visual Studio Code](https://code.visualstudio.com/) is a popular and modern choice, available for Linux, Windows and macOS.

2. Create and edit files inside the shared folder, so that they will be visible to both the host and the Docker container.

3. You will edit files from the host OS, and compile them from inside the Docker container.

---

# 2. Test the installation

1. Using your favorite text editor or IDE, move to the shared folder and create a file [`test.cpp`](test.cpp) with the following content:
   ```cpp
   #include <deal.II/base/config.h>

   #include <iostream>

   int
   main()
   {
     std::cout << "Successfully loaded deal.II." << std::endl;
     return 0;
   }
   ```
   and, in the same folder, create a file [`CMakeLists.txt`](CMakeLists.txt) with the following content:
   ```cmake
   cmake_minimum_required(VERSION 3.13.4)

   project(test_project)

   find_package(deal.II 9.5.0 REQUIRED
       HINTS $ENV{mkDealiiPrefix} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
   )
   deal_ii_initialize_cached_variables()

   add_executable(test test.cpp)
   deal_ii_setup_target(test)
   ```

2. Configure, build and run the test:
   ```bash
   mkdir build
   cd build
   cmake ..

   make -jN # N is the number of parallel processes.
   ./test
   ```
