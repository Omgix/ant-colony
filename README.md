# The Ant Colony Optimization Algorithms

This is the main source code repository for [Ant Colony Optimization Algorithms]. It contains the standard library, and documentation.

## Installing from Source
[building-from-source]: #building-from-source

### Building on *nix
1. Make sure you have installed the dependencies:

   * `g++` 4.7 or later or `clang++` 3.x or later
   * `python` 3.x (but not 2.x)
   * GNU `make` 3.81 or later
   * `cmake` 3.4.3 or later
   * `curl`
   * `git`

2. Clone the [source] with `git`:

   ```sh
   $ git clone https://github.com/Omgix/ant-colony
   $ cd ant-colony
   ```

[source]: https://github.com/Omgix/ant-colony

3. Build and install:

    ```sh
    $ cmake ./ && sudo make
    ```

    > ***Note:*** Install config can be adjusted by edit the cmake file


## Building Documentation
[building-documentation]: #building-documentation

If you’d like to build the documentation, it’s almost the same:

```sh
$ doxygen .\Doxyfile
```

The generated documentation will appear under the `ant-colony` directory for
the ABI used. 